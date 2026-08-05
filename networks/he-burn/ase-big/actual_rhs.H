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
#include <table_rates.H>
#include <temperature_table_rates.H>
#include <derived_rates.H>

using namespace amrex::literals;
using namespace ArrayUtil;

using namespace Species;
using namespace Rates;

using namespace rate_tables;

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
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 43.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Sc43) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Sc43) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Fe56) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Fe56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 44.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ti44) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ti44) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Fe55) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Fe55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mg24) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Mg24) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_C12_O16) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_C12_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 49.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cr49) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Cr49) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 11.0_rt, 22.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Na22) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Na22) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 47.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_V47) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_V47) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Co55) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Co55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_C12_C12) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_C12_C12) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 53.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Fe53) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Fe53) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 49.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cr49) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Cr49) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 13.0_rt, 27.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Al27) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Al27) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 47.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_V47) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_V47) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 12.0_rt, 25.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Mg25) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Mg25) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 22.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ne22) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ne22) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_O18) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_O18) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 52.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Mn52) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Mn52) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ni56) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ni56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_P31) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_P31) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 29.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Si29) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Si29) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 58.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ni58) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ni58) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_O16) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 52.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mn52) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Mn52) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(8.0_rt, 16.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 8.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_O16_O16) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_O16_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_P31) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_P31) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_N13) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_N13) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 29.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Si29) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Si29) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 37.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ar37) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ar37) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 19.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_F19) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_F19) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 40.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ca40) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ca40) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 38.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_K38) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_K38) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 36.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ar36) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ar36) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_F18) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_F18) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 43.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Sc43) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Sc43) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 21.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ne21) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ne21) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ne20) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ne20) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Co56) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Co56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 50.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mn50) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Mn50) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 54.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Fe54) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Fe54) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Co55) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Co55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 53.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Fe53) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Fe53) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 13.0_rt, 27.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Al27) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Al27) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 25.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mg25) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Mg25) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 35.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cl35) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Cl35) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 33.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_S33) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_S33) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 34.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cl34) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Cl34) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_C12) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_C12) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_O16) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_N14) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_N14) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 45.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ti45) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ti45) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 42.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Sc42) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Sc42) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_O17) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_O17) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_C14) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_C14) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 37.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ar37) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ar37) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_N15) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_N15) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 19.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_F19) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_F19) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_C14) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_C14) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 13.0_rt, 26.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Al26) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Al26) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 21.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ne21) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ne21) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 32.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_S32) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_S32) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 59.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cu59) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Cu59) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 30.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_P30) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_P30) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Co56) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Co56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 28.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Si28) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Si28) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 54.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Fe54) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Fe54) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 57.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ni57) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ni57) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Be8) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Be8) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 41.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ca41) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ca41) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 57.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Co57) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Co57) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 39.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_K39) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_K39) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Fe55) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Fe55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 35.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cl35) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Cl35) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_C13) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_C13) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 33.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_S33) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_S33) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_C12) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_C12) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 41.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ca41) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ca41) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 39.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_K39) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_K39) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_O18) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_O18) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 48.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cr48) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Cr48) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 45.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ti45) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ti45) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 46.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_V46) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_V46) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_O17) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_O17) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_N15) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_N15) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Mn51) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Mn51) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 23.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Na23) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Na23) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Fe52) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Fe52) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 51.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mn51) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Mn51) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_He4) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_He4) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 11.0_rt, 23.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Na23) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Na23) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_C12_Ne20) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_C12_Ne20) = dlog_scor_dT;
        }
    }

}
#endif


template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                    T& rate_eval) {

    // here Y is consistent with state.xn

    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    const tf_t tfactors = evaluate_tfactors(state.T);

    // Precompute screening terms

#ifdef SCREENING
    evaluate_screening<do_T_derivatives, T>(state, Y, rate_eval);
#endif

    // Fill in different rates

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);

    temp_tabular::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // fill modified rates next -- these can have ReacLib or
    // TemperatureTabular/StarLib rates as the original rate
    modified_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate Derived Rates next. This should go last but before
    // approx rates, since those may have ReacLibRate,
    // TemperatureTabularRate / StarLibRate, or DerivedRate as the
    // underlying rate.

    fill_derived_rates<do_T_derivatives, T>(tfactors, rate_eval);


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

    // Calculate tabular weak rates

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0_rt;

    amrex::Real log_temp = std::log10(state.T);
    amrex::Real log_rhoy = std::log10(rhoy);

    tabular_evaluate(j_Fe52_Mn52_meta, j_Fe52_Mn52_rhoy, j_Fe52_Mn52_temp, j_Fe52_Mn52_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe52_to_Mn52_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_Mn52_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Fe52) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Mn52_Fe52_meta, j_Mn52_Fe52_rhoy, j_Mn52_Fe52_temp, j_Mn52_Fe52_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Mn52_to_Fe52_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mn52_to_Fe52_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Mn52) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co55_Fe55_meta, j_Co55_Fe55_rhoy, j_Co55_Fe55_temp, j_Co55_Fe55_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co55_to_Fe55_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co55_to_Fe55_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Co55) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Fe55_Co55_meta, j_Fe55_Co55_rhoy, j_Fe55_Co55_temp, j_Fe55_Co55_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe55_to_Co55_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe55_to_Co55_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Fe55) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Fe56_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co56_to_Fe56_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Ni56_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co56_to_Ni56_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Fe56_Co56_meta, j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe56_to_Co56_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe56_to_Co56_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Fe56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ni56_Co56_meta, j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ni56_to_Co56_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_Co56_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Ni56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co57_Ni57_meta, j_Co57_Ni57_rhoy, j_Co57_Ni57_temp, j_Co57_Ni57_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co57_to_Ni57_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co57_to_Ni57_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Co57) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ni57_Co57_meta, j_Ni57_Co57_rhoy, j_Ni57_Co57_temp, j_Ni57_Co57_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ni57_to_Co57_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni57_to_Co57_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Ni57) * (edot_nu + edot_gamma);

    tabular_evaluate(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_n_to_p_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_to_p_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(N) * (edot_nu + edot_gamma);

    tabular_evaluate(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_p_to_n_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_to_n_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(H1) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F18_O18_meta, j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F18_to_O18_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_O18_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O18_F18_meta, j_O18_F18_rhoy, j_O18_F18_temp, j_O18_F18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O18_to_F18_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O18_to_F18_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(O18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Na22_Ne22_meta, j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Na22_to_Ne22_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na22_to_Ne22_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Na22) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne22_Na22_meta, j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne22_to_Na22_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne22_to_Na22_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Ne22) * (edot_nu + edot_gamma);


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

    constexpr int do_T_derivatives = 0;
    [[maybe_unused]] amrex::Real rate, drate_dT, edot_nu, edot_gamma;
    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    rate_eval.enuc_weak = 0.0_rt;

    // Compute all weak rates and get ydot_weak

    const tf_t tfactors = evaluate_tfactors(state.T);

    {
        // C14_to_N14_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_C14_to_N14_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C14_to_N14_reaclib) = rate;
    }

    {
        // N13_to_C13_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_N13_to_C13_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_N13_to_C13_reaclib) = rate;
    }

    amrex::Real log_temp = std::log10(state.T);
    amrex::Real log_rhoy = std::log10(rhoy);

    tabular_evaluate(j_Fe52_Mn52_meta, j_Fe52_Mn52_rhoy, j_Fe52_Mn52_temp, j_Fe52_Mn52_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe52_to_Mn52_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Fe52) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Mn52_Fe52_meta, j_Mn52_Fe52_rhoy, j_Mn52_Fe52_temp, j_Mn52_Fe52_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Mn52_to_Fe52_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Mn52) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co55_Fe55_meta, j_Co55_Fe55_rhoy, j_Co55_Fe55_temp, j_Co55_Fe55_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co55_to_Fe55_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Co55) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Fe55_Co55_meta, j_Fe55_Co55_rhoy, j_Fe55_Co55_temp, j_Fe55_Co55_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe55_to_Co55_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Fe55) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Fe56_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Ni56_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Fe56_Co56_meta, j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe56_to_Co56_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Fe56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ni56_Co56_meta, j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ni56_to_Co56_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Ni56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co57_Ni57_meta, j_Co57_Ni57_rhoy, j_Co57_Ni57_temp, j_Co57_Ni57_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co57_to_Ni57_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Co57) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ni57_Co57_meta, j_Ni57_Co57_rhoy, j_Ni57_Co57_temp, j_Ni57_Co57_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ni57_to_Co57_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Ni57) * (edot_nu + edot_gamma);

    tabular_evaluate(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_n_to_p_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(N) * (edot_nu + edot_gamma);

    tabular_evaluate(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_p_to_n_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(H1) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F18_O18_meta, j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F18_to_O18_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O18_F18_meta, j_O18_F18_rhoy, j_O18_F18_temp, j_O18_F18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O18_to_F18_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(O18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Na22_Ne22_meta, j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Na22_to_Ne22_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Na22) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne22_Na22_meta, j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dT, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne22_to_Na22_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Ne22) * (edot_nu + edot_gamma);


    const auto& screened_rates = rate_eval.screened_rates;

    ydot_nuc(N) =
        (-screened_rates(k_n_to_p_weaktab)*Y(N) + screened_rates(k_p_to_n_weaktab)*Y(H1));

    ydot_nuc(H1) =
        (screened_rates(k_n_to_p_weaktab)*Y(N) + -screened_rates(k_p_to_n_weaktab)*Y(H1));

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(C12) = 0.0_rt;

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_reaclib)*Y(N13);

    ydot_nuc(C14) =
        -screened_rates(k_C14_to_N14_reaclib)*Y(C14);

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_reaclib)*Y(N13);

    ydot_nuc(N14) =
        screened_rates(k_C14_to_N14_reaclib)*Y(C14);

    ydot_nuc(N15) = 0.0_rt;

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(O17) = 0.0_rt;

    ydot_nuc(O18) =
        (screened_rates(k_F18_to_O18_weaktab)*Y(F18) + -screened_rates(k_O18_to_F18_weaktab)*Y(O18));

    ydot_nuc(F18) =
        (-screened_rates(k_F18_to_O18_weaktab)*Y(F18) + screened_rates(k_O18_to_F18_weaktab)*Y(O18));

    ydot_nuc(F19) = 0.0_rt;

    ydot_nuc(Ne20) = 0.0_rt;

    ydot_nuc(Ne21) = 0.0_rt;

    ydot_nuc(Ne22) =
        (screened_rates(k_Na22_to_Ne22_weaktab)*Y(Na22) + -screened_rates(k_Ne22_to_Na22_weaktab)*Y(Ne22));

    ydot_nuc(Na22) =
        (-screened_rates(k_Na22_to_Ne22_weaktab)*Y(Na22) + screened_rates(k_Ne22_to_Na22_weaktab)*Y(Ne22));

    ydot_nuc(Na23) = 0.0_rt;

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Mg25) = 0.0_rt;

    ydot_nuc(Al26) = 0.0_rt;

    ydot_nuc(Al27) = 0.0_rt;

    ydot_nuc(Si28) = 0.0_rt;

    ydot_nuc(Si29) = 0.0_rt;

    ydot_nuc(P30) = 0.0_rt;

    ydot_nuc(P31) = 0.0_rt;

    ydot_nuc(S32) = 0.0_rt;

    ydot_nuc(S33) = 0.0_rt;

    ydot_nuc(Cl34) = 0.0_rt;

    ydot_nuc(Cl35) = 0.0_rt;

    ydot_nuc(Ar36) = 0.0_rt;

    ydot_nuc(Ar37) = 0.0_rt;

    ydot_nuc(K38) = 0.0_rt;

    ydot_nuc(K39) = 0.0_rt;

    ydot_nuc(Ca40) = 0.0_rt;

    ydot_nuc(Ca41) = 0.0_rt;

    ydot_nuc(Sc42) = 0.0_rt;

    ydot_nuc(Sc43) = 0.0_rt;

    ydot_nuc(Ti44) = 0.0_rt;

    ydot_nuc(Ti45) = 0.0_rt;

    ydot_nuc(V46) = 0.0_rt;

    ydot_nuc(V47) = 0.0_rt;

    ydot_nuc(Cr48) = 0.0_rt;

    ydot_nuc(Cr49) = 0.0_rt;

    ydot_nuc(Mn50) = 0.0_rt;

    ydot_nuc(Mn51) = 0.0_rt;

    ydot_nuc(Mn52) =
        (screened_rates(k_Fe52_to_Mn52_weaktab)*Y(Fe52) + -screened_rates(k_Mn52_to_Fe52_weaktab)*Y(Mn52));

    ydot_nuc(Fe52) =
        (-screened_rates(k_Fe52_to_Mn52_weaktab)*Y(Fe52) + screened_rates(k_Mn52_to_Fe52_weaktab)*Y(Mn52));

    ydot_nuc(Fe53) = 0.0_rt;

    ydot_nuc(Fe54) = 0.0_rt;

    ydot_nuc(Fe55) =
        (screened_rates(k_Co55_to_Fe55_weaktab)*Y(Co55) + -screened_rates(k_Fe55_to_Co55_weaktab)*Y(Fe55));

    ydot_nuc(Fe56) =
        (screened_rates(k_Co56_to_Fe56_weaktab)*Y(Co56) + -screened_rates(k_Fe56_to_Co56_weaktab)*Y(Fe56));

    ydot_nuc(Co55) =
        (-screened_rates(k_Co55_to_Fe55_weaktab)*Y(Co55) + screened_rates(k_Fe55_to_Co55_weaktab)*Y(Fe55));

    ydot_nuc(Co56) =
        (-screened_rates(k_Co56_to_Fe56_weaktab)*Y(Co56) + screened_rates(k_Fe56_to_Co56_weaktab)*Y(Fe56)) +
        (screened_rates(k_Ni56_to_Co56_weaktab)*Y(Ni56) + -screened_rates(k_Co56_to_Ni56_weaktab)*Y(Co56));

    ydot_nuc(Co57) =
        (screened_rates(k_Ni57_to_Co57_weaktab)*Y(Ni57) + -screened_rates(k_Co57_to_Ni57_weaktab)*Y(Co57));

    ydot_nuc(Ni56) =
        (-screened_rates(k_Ni56_to_Co56_weaktab)*Y(Ni56) + screened_rates(k_Co56_to_Ni56_weaktab)*Y(Co56));

    ydot_nuc(Ni57) =
        (-screened_rates(k_Ni57_to_Co57_weaktab)*Y(Ni57) + screened_rates(k_Co57_to_Ni57_weaktab)*Y(Co57));

    ydot_nuc(Ni58) = 0.0_rt;

    ydot_nuc(Cu59) = 0.0_rt;

    ydot_nuc(Zn60) = 0.0_rt;

    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        (-screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(Mn50)*Y(N)*state.rho + screened_rates(k_Mn51_to_n_Mn50_derived)*Y(Mn51)) +
        (-screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(Mn51)*Y(N)*state.rho + screened_rates(k_Mn52_to_n_Mn51_derived)*Y(Mn52)) +
        (-screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(Fe52)*Y(N)*state.rho + screened_rates(k_Fe53_to_n_Fe52_derived)*Y(Fe53)) +
        (-screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(Fe53)*Y(N)*state.rho + screened_rates(k_Fe54_to_n_Fe53_derived)*Y(Fe54)) +
        (-screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(Fe54)*Y(N)*state.rho + screened_rates(k_Fe55_to_n_Fe54_derived)*Y(Fe55)) +
        (-screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(Fe55)*Y(N)*state.rho + screened_rates(k_Fe56_to_n_Fe55_derived)*Y(Fe56)) +
        (-screened_rates(k_n_Co55_to_Co56_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_Co56_to_n_Co55_derived)*Y(Co56)) +
        (-screened_rates(k_n_Co56_to_Co57_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Co57_to_n_Co56_derived)*Y(Co57)) +
        (-screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_Ni57_to_n_Ni56_derived)*Y(Ni57)) +
        (-screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(Ni57)*Y(N)*state.rho + screened_rates(k_Ni58_to_n_Ni57_derived)*Y(Ni58)) +
        (-screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*Y(N)*state.rho + screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*Y(Mn52)*state.rho) +
        (-screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*Y(He4)*state.rho) +
        (-screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*Y(N)*state.rho + screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*Y(N)*state.rho + screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (-screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*Y(N)*state.rho + screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*Y(He4)*state.rho) +
        (-screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*Y(N)*state.rho + screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*Y(Ni57)*state.rho) +
        (-screened_rates(k_n_to_p_weaktab)*Y(N) + screened_rates(k_p_to_n_weaktab)*Y(H1));

    ydot_nuc(H1) =
        (-screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13_derived)*Y(N14)) +
        (-screened_rates(k_p_C14_to_N15_reaclib)*Y(C14)*Y(H1)*state.rho + screened_rates(k_N15_to_p_C14_derived)*Y(N15)) +
        (-screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (-screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(Ne22)*Y(H1)*state.rho + screened_rates(k_Na23_to_p_Ne22_derived)*Y(Na23)) +
        (-screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        (-screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(Mg25)*Y(H1)*state.rho + screened_rates(k_Al26_to_p_Mg25_derived)*Y(Al26)) +
        (-screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        (-screened_rates(k_p_Si29_to_P30_reaclib)*Y(Si29)*Y(H1)*state.rho + screened_rates(k_P30_to_p_Si29_derived)*Y(P30)) +
        (-screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        (-screened_rates(k_p_S33_to_Cl34_reaclib)*Y(S33)*Y(H1)*state.rho + screened_rates(k_Cl34_to_p_S33_derived)*Y(Cl34)) +
        (-screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_Ar36_to_p_Cl35_derived)*Y(Ar36)) +
        (-screened_rates(k_p_Ar37_to_K38_reaclib)*Y(Ar37)*Y(H1)*state.rho + screened_rates(k_K38_to_p_Ar37_derived)*Y(K38)) +
        (-screened_rates(k_p_K39_to_Ca40_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ca40_to_p_K39_derived)*Y(Ca40)) +
        (-screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(Ca41)*Y(H1)*state.rho + screened_rates(k_Sc42_to_p_Ca41_derived)*Y(Sc42)) +
        (-screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ti44_to_p_Sc43_derived)*Y(Ti44)) +
        (-screened_rates(k_p_Ti45_to_V46_reaclib)*Y(Ti45)*Y(H1)*state.rho + screened_rates(k_V46_to_p_Ti45_derived)*Y(V46)) +
        (-screened_rates(k_p_V47_to_Cr48_reaclib)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Cr48_to_p_V47_derived)*Y(Cr48)) +
        (-screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(Cr49)*Y(H1)*state.rho + screened_rates(k_Mn50_to_p_Cr49_derived)*Y(Mn50)) +
        (-screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_Fe52_to_p_Mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55_derived)*Y(Ni56)) +
        (0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_p_Na23_to_C12_C12_derived)*Y(Na23)*Y(H1)*state.rho) +
        (screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_p_P31_to_O16_O16_derived)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (-screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*Y(He4)*state.rho + -screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*Y(H1)*state.rho) +
        (screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*Y(Na22)*state.rho + -screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*Y(He4)*state.rho + -screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*Y(P30)*state.rho + -screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*Y(H1)*state.rho) +
        (-screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*Y(He4)*state.rho + -screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*Y(S32)*state.rho) +
        (screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*Y(K38)*state.rho + -screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*Y(H1)*state.rho) +
        (-screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*Y(He4)*state.rho) +
        (screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*Y(Sc42)*state.rho + -screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*Y(He4)*state.rho) +
        (screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*Y(V46)*state.rho + -screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*Y(He4)*state.rho + -screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(Mn52)*Y(H1)*state.rho + screened_rates(k_Fe53_to_p_Mn52_derived)*Y(Fe53)) +
        (-screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(Fe54)*Y(H1)*state.rho + screened_rates(k_Co55_to_p_Fe54_derived)*Y(Co55)) +
        (-screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(Fe55)*Y(H1)*state.rho + screened_rates(k_Co56_to_p_Fe55_derived)*Y(Co56)) +
        (-screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(Fe56)*Y(H1)*state.rho + screened_rates(k_Co57_to_p_Fe56_derived)*Y(Co57)) +
        (-screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(Co56)*Y(H1)*state.rho + screened_rates(k_Ni57_to_p_Co56_derived)*Y(Ni57)) +
        (-screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_Ni58_to_p_Co57_derived)*Y(Ni58)) +
        (-screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(Ni58)*Y(H1)*state.rho + screened_rates(k_Cu59_to_p_Ni58_derived)*Y(Cu59)) +
        (-screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_Zn60_to_p_Cu59_derived)*Y(Zn60)) +
        (screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*Y(Mn50)*state.rho + -screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*Y(Mn52)*state.rho + -screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*Y(N)*state.rho + -screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*Y(He4)*state.rho + -screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*Y(N)*state.rho + -screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*Y(Ni56)*state.rho) +
        (screened_rates(k_n_to_p_weaktab)*Y(N) + -screened_rates(k_p_to_n_weaktab)*Y(H1)) +
        (-screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*Y(H1)*state.rho + screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_He4_C14_to_O18_reaclib)*Y(C14)*Y(He4)*state.rho + screened_rates(k_O18_to_He4_C14_derived)*Y(O18)) +
        (-screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (-screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (-screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*Y(O17)*state.rho + screened_rates(k_Ne21_to_He4_O17_derived)*Y(Ne21)) +
        (-screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(He4)*Y(O18)*state.rho + screened_rates(k_Ne22_to_He4_O18_derived)*Y(Ne22)) +
        (-screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (-screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*Y(He4)*state.rho + screened_rates(k_Na23_to_He4_F19_derived)*Y(Na23)) +
        (-screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (-screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(He4)*Y(Ne21)*state.rho + screened_rates(k_Mg25_to_He4_Ne21_derived)*Y(Mg25)) +
        (-screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(He4)*Y(Na22)*state.rho + screened_rates(k_Al26_to_He4_Na22_derived)*Y(Al26)) +
        (-screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(He4)*Y(Na23)*state.rho + screened_rates(k_Al27_to_He4_Na23_derived)*Y(Al27)) +
        (-screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        (-screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(He4)*Y(Mg25)*state.rho + screened_rates(k_Si29_to_He4_Mg25_derived)*Y(Si29)) +
        (-screened_rates(k_He4_Al26_to_P30_reaclib)*Y(Al26)*Y(He4)*state.rho + screened_rates(k_P30_to_He4_Al26_derived)*Y(P30)) +
        (-screened_rates(k_He4_Al27_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27_derived)*Y(P31)) +
        (-screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        (-screened_rates(k_He4_Si29_to_S33_reaclib)*Y(He4)*Y(Si29)*state.rho + screened_rates(k_S33_to_He4_Si29_derived)*Y(S33)) +
        (-screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(He4)*Y(P30)*state.rho + screened_rates(k_Cl34_to_He4_P30_derived)*Y(Cl34)) +
        (-screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(He4)*Y(P31)*state.rho + screened_rates(k_Cl35_to_He4_P31_derived)*Y(Cl35)) +
        (-screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_He4_S32_derived)*Y(Ar36)) +
        (-screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(He4)*Y(S33)*state.rho + screened_rates(k_Ar37_to_He4_S33_derived)*Y(Ar37)) +
        (-screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(Cl34)*Y(He4)*state.rho + screened_rates(k_K38_to_He4_Cl34_derived)*Y(K38)) +
        (-screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35_derived)*Y(K39)) +
        (-screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_He4_Ar36_derived)*Y(Ca40)) +
        (-screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(Ar37)*Y(He4)*state.rho + screened_rates(k_Ca41_to_He4_Ar37_derived)*Y(Ca41)) +
        (-screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(He4)*Y(K38)*state.rho + screened_rates(k_Sc42_to_He4_K38_derived)*Y(Sc42)) +
        (-screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(He4)*Y(K39)*state.rho + screened_rates(k_Sc43_to_He4_K39_derived)*Y(Sc43)) +
        (-screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40_derived)*Y(Ti44)) +
        (-screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(Ca41)*Y(He4)*state.rho + screened_rates(k_Ti45_to_He4_Ca41_derived)*Y(Ti45)) +
        (-screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(He4)*Y(Sc42)*state.rho + screened_rates(k_V46_to_He4_Sc42_derived)*Y(V46)) +
        (-screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_V47_to_He4_Sc43_derived)*Y(V47)) +
        (-screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_He4_Ti44_derived)*Y(Cr48)) +
        (-screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(He4)*Y(Ti45)*state.rho + screened_rates(k_Cr49_to_He4_Ti45_derived)*Y(Cr49)) +
        (-screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(He4)*Y(V46)*state.rho + screened_rates(k_Mn50_to_He4_V46_derived)*Y(Mn50)) +
        (-screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Mn51_to_He4_V47_derived)*Y(Mn51)) +
        (-screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51_derived)*Y(Co55)) +
        (-screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (-screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (-screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*Y(He4)*state.rho + screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*Y(H1)*state.rho) +
        (screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*Y(Na22)*state.rho + screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*Y(H1)*state.rho) +
        (screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*Y(He4)*state.rho + screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*Y(H1)*state.rho) +
        (screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*Y(P30)*state.rho + screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*Y(H1)*state.rho) +
        (screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*Y(He4)*state.rho + screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*Y(H1)*state.rho) +
        (screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*Y(S32)*state.rho) +
        (-screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*Y(K38)*state.rho + screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*Y(H1)*state.rho) +
        (screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*Y(He4)*state.rho) +
        (-screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*Y(Sc42)*state.rho + screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*Y(H1)*state.rho) +
        (screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*Y(He4)*state.rho) +
        (-screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*Y(V46)*state.rho + screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*Y(He4)*state.rho + screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + 3.0*screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12)) +
        (-screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(He4)*Y(Mn52)*state.rho + screened_rates(k_Co56_to_He4_Mn52_derived)*Y(Co56)) +
        (-screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_Ni57_to_He4_Fe53_derived)*Y(Ni57)) +
        (-screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_Ni58_to_He4_Fe54_derived)*Y(Ni58)) +
        (-screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_Cu59_to_He4_Co55_derived)*Y(Cu59)) +
        (-screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(He4)*Y(Ni56)*state.rho + screened_rates(k_Zn60_to_He4_Ni56_derived)*Y(Zn60)) +
        (-screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*Y(Mn50)*state.rho + screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*Y(Mn52)*state.rho + screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*Y(Mn52)*state.rho) +
        (-screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*Y(He4)*state.rho) +
        (screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*Y(N)*state.rho + -screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*Y(N)*state.rho + -screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*Y(He4)*state.rho) +
        (screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + -screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*Y(He4)*state.rho) +
        (screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*Y(Ni56)*state.rho) +
        (screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*Y(N)*state.rho + -screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*Y(Ni57)*state.rho) +
        (screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho);

    ydot_nuc(C12) =
        (-screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + 2.0*screened_rates(k_p_Na23_to_C12_C12_derived)*Y(Na23)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (-screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*Y(Si28)*state.rho) +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12));

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_reaclib)*Y(N13) +
        (-screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13_derived)*Y(N14));

    ydot_nuc(C14) =
        -screened_rates(k_C14_to_N14_reaclib)*Y(C14) +
        (-screened_rates(k_p_C14_to_N15_reaclib)*Y(C14)*Y(H1)*state.rho + screened_rates(k_N15_to_p_C14_derived)*Y(N15)) +
        (-screened_rates(k_He4_C14_to_O18_reaclib)*Y(C14)*Y(He4)*state.rho + screened_rates(k_O18_to_He4_C14_derived)*Y(O18));

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_reaclib)*Y(N13) +
        (screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + -screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho);

    ydot_nuc(N14) =
        screened_rates(k_C14_to_N14_reaclib)*Y(C14) +
        (screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + -screened_rates(k_N14_to_p_C13_derived)*Y(N14)) +
        (-screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho);

    ydot_nuc(N15) =
        (screened_rates(k_p_C14_to_N15_reaclib)*Y(C14)*Y(H1)*state.rho + -screened_rates(k_N15_to_p_C14_derived)*Y(N15)) +
        (-screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (-screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho);

    ydot_nuc(O16) =
        (screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + 2.0*screened_rates(k_p_P31_to_O16_O16_derived)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + 2.0*screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho);

    ydot_nuc(O17) =
        (-screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*Y(O17)*state.rho + screened_rates(k_Ne21_to_He4_O17_derived)*Y(Ne21)) +
        (-screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho);

    ydot_nuc(O18) =
        (screened_rates(k_He4_C14_to_O18_reaclib)*Y(C14)*Y(He4)*state.rho + -screened_rates(k_O18_to_He4_C14_derived)*Y(O18)) +
        (-screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(He4)*Y(O18)*state.rho + screened_rates(k_Ne22_to_He4_O18_derived)*Y(Ne22)) +
        (-screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_F18_to_O18_weaktab)*Y(F18) + -screened_rates(k_O18_to_F18_weaktab)*Y(O18));

    ydot_nuc(F18) =
        (screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (-screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_F18_to_O18_weaktab)*Y(F18) + screened_rates(k_O18_to_F18_weaktab)*Y(O18));

    ydot_nuc(F19) =
        (screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + -screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*Y(He4)*state.rho + screened_rates(k_Na23_to_He4_F19_derived)*Y(Na23)) +
        (-screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (-screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*Y(He4)*state.rho + screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*Y(H1)*state.rho);

    ydot_nuc(Ne20) =
        (screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Ne21) =
        (screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*Y(O17)*state.rho + -screened_rates(k_Ne21_to_He4_O17_derived)*Y(Ne21)) +
        (-screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (-screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(He4)*Y(Ne21)*state.rho + screened_rates(k_Mg25_to_He4_Ne21_derived)*Y(Mg25)) +
        (screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Ne22) =
        (screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(He4)*Y(O18)*state.rho + -screened_rates(k_Ne22_to_He4_O18_derived)*Y(Ne22)) +
        (-screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(Ne22)*Y(H1)*state.rho + screened_rates(k_Na23_to_p_Ne22_derived)*Y(Na23)) +
        (screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*Y(He4)*state.rho + -screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*Y(H1)*state.rho) +
        (screened_rates(k_Na22_to_Ne22_weaktab)*Y(Na22) + -screened_rates(k_Ne22_to_Na22_weaktab)*Y(Ne22));

    ydot_nuc(Na22) =
        (screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*Y(H1)*state.rho + -screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (-screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(He4)*Y(Na22)*state.rho + screened_rates(k_Al26_to_He4_Na22_derived)*Y(Al26)) +
        (-screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*Y(Na22)*state.rho + screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*Y(H1)*state.rho) +
        (-screened_rates(k_Na22_to_Ne22_weaktab)*Y(Na22) + screened_rates(k_Ne22_to_Na22_weaktab)*Y(Ne22));

    ydot_nuc(Na23) =
        (screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*Y(He4)*state.rho + -screened_rates(k_Na23_to_He4_F19_derived)*Y(Na23)) +
        (screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(Ne22)*Y(H1)*state.rho + -screened_rates(k_Na23_to_p_Ne22_derived)*Y(Na23)) +
        (-screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        (-screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(He4)*Y(Na23)*state.rho + screened_rates(k_Al27_to_He4_Na23_derived)*Y(Al27)) +
        (0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_p_Na23_to_C12_C12_derived)*Y(Na23)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Mg24) =
        (screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        (-screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        (screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Mg25) =
        (screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(He4)*Y(Ne21)*state.rho + -screened_rates(k_Mg25_to_He4_Ne21_derived)*Y(Mg25)) +
        (-screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(Mg25)*Y(H1)*state.rho + screened_rates(k_Al26_to_p_Mg25_derived)*Y(Al26)) +
        (-screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(He4)*Y(Mg25)*state.rho + screened_rates(k_Si29_to_He4_Mg25_derived)*Y(Si29)) +
        (screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*Y(Na22)*state.rho + -screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*Y(H1)*state.rho);

    ydot_nuc(Al26) =
        (screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(He4)*Y(Na22)*state.rho + -screened_rates(k_Al26_to_He4_Na22_derived)*Y(Al26)) +
        (screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(Mg25)*Y(H1)*state.rho + -screened_rates(k_Al26_to_p_Mg25_derived)*Y(Al26)) +
        (-screened_rates(k_He4_Al26_to_P30_reaclib)*Y(Al26)*Y(He4)*state.rho + screened_rates(k_P30_to_He4_Al26_derived)*Y(P30)) +
        (-screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*Y(He4)*state.rho + screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*Y(H1)*state.rho);

    ydot_nuc(Al27) =
        (screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(He4)*Y(Na23)*state.rho + -screened_rates(k_Al27_to_He4_Na23_derived)*Y(Al27)) +
        (-screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        (-screened_rates(k_He4_Al27_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27_derived)*Y(P31)) +
        (screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Si28) =
        (screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + -screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        (screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        (-screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        (0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(Si29) =
        (screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(He4)*Y(Mg25)*state.rho + -screened_rates(k_Si29_to_He4_Mg25_derived)*Y(Si29)) +
        (-screened_rates(k_p_Si29_to_P30_reaclib)*Y(Si29)*Y(H1)*state.rho + screened_rates(k_P30_to_p_Si29_derived)*Y(P30)) +
        (-screened_rates(k_He4_Si29_to_S33_reaclib)*Y(He4)*Y(Si29)*state.rho + screened_rates(k_S33_to_He4_Si29_derived)*Y(S33)) +
        (screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*Y(He4)*state.rho + -screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*Y(H1)*state.rho);

    ydot_nuc(P30) =
        (screened_rates(k_He4_Al26_to_P30_reaclib)*Y(Al26)*Y(He4)*state.rho + -screened_rates(k_P30_to_He4_Al26_derived)*Y(P30)) +
        (screened_rates(k_p_Si29_to_P30_reaclib)*Y(Si29)*Y(H1)*state.rho + -screened_rates(k_P30_to_p_Si29_derived)*Y(P30)) +
        (-screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(He4)*Y(P30)*state.rho + screened_rates(k_Cl34_to_He4_P30_derived)*Y(Cl34)) +
        (-screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*Y(P30)*state.rho + screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*Y(H1)*state.rho);

    ydot_nuc(P31) =
        (screened_rates(k_He4_Al27_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho + -screened_rates(k_P31_to_He4_Al27_derived)*Y(P31)) +
        (-screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        (-screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(He4)*Y(P31)*state.rho + screened_rates(k_Cl35_to_He4_P31_derived)*Y(Cl35)) +
        (0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_p_P31_to_O16_O16_derived)*Y(P31)*Y(H1)*state.rho) +
        (screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(S32) =
        (screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + -screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        (screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        (-screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_He4_S32_derived)*Y(Ar36)) +
        (screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*Y(S32)*state.rho);

    ydot_nuc(S33) =
        (screened_rates(k_He4_Si29_to_S33_reaclib)*Y(He4)*Y(Si29)*state.rho + -screened_rates(k_S33_to_He4_Si29_derived)*Y(S33)) +
        (-screened_rates(k_p_S33_to_Cl34_reaclib)*Y(S33)*Y(H1)*state.rho + screened_rates(k_Cl34_to_p_S33_derived)*Y(Cl34)) +
        (-screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(He4)*Y(S33)*state.rho + screened_rates(k_Ar37_to_He4_S33_derived)*Y(Ar37)) +
        (screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*Y(P30)*state.rho + -screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*Y(H1)*state.rho);

    ydot_nuc(Cl34) =
        (screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(He4)*Y(P30)*state.rho + -screened_rates(k_Cl34_to_He4_P30_derived)*Y(Cl34)) +
        (screened_rates(k_p_S33_to_Cl34_reaclib)*Y(S33)*Y(H1)*state.rho + -screened_rates(k_Cl34_to_p_S33_derived)*Y(Cl34)) +
        (-screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(Cl34)*Y(He4)*state.rho + screened_rates(k_K38_to_He4_Cl34_derived)*Y(K38)) +
        (-screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*Y(He4)*state.rho + screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*Y(H1)*state.rho);

    ydot_nuc(Cl35) =
        (screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(He4)*Y(P31)*state.rho + -screened_rates(k_Cl35_to_He4_P31_derived)*Y(Cl35)) +
        (-screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_Ar36_to_p_Cl35_derived)*Y(Ar36)) +
        (-screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35_derived)*Y(K39)) +
        (-screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*Y(S32)*state.rho);

    ydot_nuc(Ar36) =
        (screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(He4)*Y(S32)*state.rho + -screened_rates(k_Ar36_to_He4_S32_derived)*Y(Ar36)) +
        (screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_Ar36_to_p_Cl35_derived)*Y(Ar36)) +
        (-screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_He4_Ar36_derived)*Y(Ca40)) +
        (screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*Y(He4)*state.rho);

    ydot_nuc(Ar37) =
        (screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(He4)*Y(S33)*state.rho + -screened_rates(k_Ar37_to_He4_S33_derived)*Y(Ar37)) +
        (-screened_rates(k_p_Ar37_to_K38_reaclib)*Y(Ar37)*Y(H1)*state.rho + screened_rates(k_K38_to_p_Ar37_derived)*Y(K38)) +
        (-screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(Ar37)*Y(He4)*state.rho + screened_rates(k_Ca41_to_He4_Ar37_derived)*Y(Ca41)) +
        (screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*Y(He4)*state.rho + -screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*Y(H1)*state.rho);

    ydot_nuc(K38) =
        (screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(Cl34)*Y(He4)*state.rho + -screened_rates(k_K38_to_He4_Cl34_derived)*Y(K38)) +
        (screened_rates(k_p_Ar37_to_K38_reaclib)*Y(Ar37)*Y(H1)*state.rho + -screened_rates(k_K38_to_p_Ar37_derived)*Y(K38)) +
        (-screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(He4)*Y(K38)*state.rho + screened_rates(k_Sc42_to_He4_K38_derived)*Y(Sc42)) +
        (-screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*Y(K38)*state.rho + screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*Y(H1)*state.rho);

    ydot_nuc(K39) =
        (screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(Cl35)*Y(He4)*state.rho + -screened_rates(k_K39_to_He4_Cl35_derived)*Y(K39)) +
        (-screened_rates(k_p_K39_to_Ca40_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ca40_to_p_K39_derived)*Y(Ca40)) +
        (-screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(He4)*Y(K39)*state.rho + screened_rates(k_Sc43_to_He4_K39_derived)*Y(Sc43)) +
        (-screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*Y(He4)*state.rho);

    ydot_nuc(Ca40) =
        (screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(Ar36)*Y(He4)*state.rho + -screened_rates(k_Ca40_to_He4_Ar36_derived)*Y(Ca40)) +
        (screened_rates(k_p_K39_to_Ca40_reaclib)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_Ca40_to_p_K39_derived)*Y(Ca40)) +
        (-screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40_derived)*Y(Ti44)) +
        (screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*Y(He4)*state.rho);

    ydot_nuc(Ca41) =
        (screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(Ar37)*Y(He4)*state.rho + -screened_rates(k_Ca41_to_He4_Ar37_derived)*Y(Ca41)) +
        (-screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(Ca41)*Y(H1)*state.rho + screened_rates(k_Sc42_to_p_Ca41_derived)*Y(Sc42)) +
        (-screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(Ca41)*Y(He4)*state.rho + screened_rates(k_Ti45_to_He4_Ca41_derived)*Y(Ti45)) +
        (screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*Y(K38)*state.rho + -screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*Y(H1)*state.rho);

    ydot_nuc(Sc42) =
        (screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(He4)*Y(K38)*state.rho + -screened_rates(k_Sc42_to_He4_K38_derived)*Y(Sc42)) +
        (screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(Ca41)*Y(H1)*state.rho + -screened_rates(k_Sc42_to_p_Ca41_derived)*Y(Sc42)) +
        (-screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(He4)*Y(Sc42)*state.rho + screened_rates(k_V46_to_He4_Sc42_derived)*Y(V46)) +
        (-screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*Y(Sc42)*state.rho + screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*Y(H1)*state.rho);

    ydot_nuc(Sc43) =
        (screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(He4)*Y(K39)*state.rho + -screened_rates(k_Sc43_to_He4_K39_derived)*Y(Sc43)) +
        (-screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ti44_to_p_Sc43_derived)*Y(Ti44)) +
        (-screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_V47_to_He4_Sc43_derived)*Y(V47)) +
        (-screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*Y(He4)*state.rho);

    ydot_nuc(Ti44) =
        (screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(Ca40)*Y(He4)*state.rho + -screened_rates(k_Ti44_to_He4_Ca40_derived)*Y(Ti44)) +
        (screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_Ti44_to_p_Sc43_derived)*Y(Ti44)) +
        (-screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_He4_Ti44_derived)*Y(Cr48)) +
        (screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho);

    ydot_nuc(Ti45) =
        (screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(Ca41)*Y(He4)*state.rho + -screened_rates(k_Ti45_to_He4_Ca41_derived)*Y(Ti45)) +
        (-screened_rates(k_p_Ti45_to_V46_reaclib)*Y(Ti45)*Y(H1)*state.rho + screened_rates(k_V46_to_p_Ti45_derived)*Y(V46)) +
        (-screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(He4)*Y(Ti45)*state.rho + screened_rates(k_Cr49_to_He4_Ti45_derived)*Y(Cr49)) +
        (screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*Y(Sc42)*state.rho + -screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*Y(H1)*state.rho);

    ydot_nuc(V46) =
        (screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(He4)*Y(Sc42)*state.rho + -screened_rates(k_V46_to_He4_Sc42_derived)*Y(V46)) +
        (screened_rates(k_p_Ti45_to_V46_reaclib)*Y(Ti45)*Y(H1)*state.rho + -screened_rates(k_V46_to_p_Ti45_derived)*Y(V46)) +
        (-screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(He4)*Y(V46)*state.rho + screened_rates(k_Mn50_to_He4_V46_derived)*Y(Mn50)) +
        (-screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*Y(V46)*state.rho + screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*Y(H1)*state.rho);

    ydot_nuc(V47) =
        (screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(He4)*Y(Sc43)*state.rho + -screened_rates(k_V47_to_He4_Sc43_derived)*Y(V47)) +
        (-screened_rates(k_p_V47_to_Cr48_reaclib)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Cr48_to_p_V47_derived)*Y(Cr48)) +
        (-screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Mn51_to_He4_V47_derived)*Y(Mn51)) +
        (-screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*Y(H1)*state.rho + screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho);

    ydot_nuc(Cr48) =
        (screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_Cr48_to_He4_Ti44_derived)*Y(Cr48)) +
        (screened_rates(k_p_V47_to_Cr48_reaclib)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_Cr48_to_p_V47_derived)*Y(Cr48)) +
        (-screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(H1)*state.rho);

    ydot_nuc(Cr49) =
        (screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(He4)*Y(Ti45)*state.rho + -screened_rates(k_Cr49_to_He4_Ti45_derived)*Y(Cr49)) +
        (-screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(Cr49)*Y(H1)*state.rho + screened_rates(k_Mn50_to_p_Cr49_derived)*Y(Mn50)) +
        (screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*Y(V46)*state.rho + -screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*Y(He4)*state.rho + screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*Y(H1)*state.rho);

    ydot_nuc(Mn50) =
        (screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(He4)*Y(V46)*state.rho + -screened_rates(k_Mn50_to_He4_V46_derived)*Y(Mn50)) +
        (screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(Cr49)*Y(H1)*state.rho + -screened_rates(k_Mn50_to_p_Cr49_derived)*Y(Mn50)) +
        (-screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(Mn50)*Y(N)*state.rho + screened_rates(k_Mn51_to_n_Mn50_derived)*Y(Mn51)) +
        (-screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*Y(Mn50)*state.rho + screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*Y(H1)*state.rho);

    ydot_nuc(Mn51) =
        (screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(He4)*Y(V47)*state.rho + -screened_rates(k_Mn51_to_He4_V47_derived)*Y(Mn51)) +
        (-screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_Fe52_to_p_Mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51_derived)*Y(Co55)) +
        (screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(H1)*state.rho) +
        (screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(Mn50)*Y(N)*state.rho + -screened_rates(k_Mn51_to_n_Mn50_derived)*Y(Mn51)) +
        (-screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(Mn51)*Y(N)*state.rho + screened_rates(k_Mn52_to_n_Mn51_derived)*Y(Mn52)) +
        (-screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(H1)*state.rho);

    ydot_nuc(Mn52) =
        (screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*Y(He4)*state.rho + -screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*Y(H1)*state.rho) +
        (screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(Mn51)*Y(N)*state.rho + -screened_rates(k_Mn52_to_n_Mn51_derived)*Y(Mn52)) +
        (-screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(Mn52)*Y(H1)*state.rho + screened_rates(k_Fe53_to_p_Mn52_derived)*Y(Fe53)) +
        (-screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(He4)*Y(Mn52)*state.rho + screened_rates(k_Co56_to_He4_Mn52_derived)*Y(Co56)) +
        (-screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*Y(Mn52)*state.rho + screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*Y(N)*state.rho + -screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*Y(H1)*state.rho) +
        (screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*Y(Mn52)*state.rho) +
        (screened_rates(k_Fe52_to_Mn52_weaktab)*Y(Fe52) + -screened_rates(k_Mn52_to_Fe52_weaktab)*Y(Mn52));

    ydot_nuc(Fe52) =
        (screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Fe52_to_He4_Cr48_derived)*Y(Fe52)) +
        (screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(Mn51)*Y(H1)*state.rho + -screened_rates(k_Fe52_to_p_Mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (-screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(Fe52)*Y(N)*state.rho + screened_rates(k_Fe53_to_n_Fe52_derived)*Y(Fe53)) +
        (-screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*Y(N)*state.rho + screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*Y(H1)*state.rho) +
        (-screened_rates(k_Fe52_to_Mn52_weaktab)*Y(Fe52) + screened_rates(k_Mn52_to_Fe52_weaktab)*Y(Mn52));

    ydot_nuc(Fe53) =
        (screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(Mn52)*Y(H1)*state.rho + -screened_rates(k_Fe53_to_p_Mn52_derived)*Y(Fe53)) +
        (screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(Fe52)*Y(N)*state.rho + -screened_rates(k_Fe53_to_n_Fe52_derived)*Y(Fe53)) +
        (-screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(Fe53)*Y(N)*state.rho + screened_rates(k_Fe54_to_n_Fe53_derived)*Y(Fe54)) +
        (-screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_Ni57_to_He4_Fe53_derived)*Y(Ni57)) +
        (screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*Y(Mn50)*state.rho + -screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*Y(He4)*state.rho);

    ydot_nuc(Fe54) =
        (screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(Fe53)*Y(N)*state.rho + -screened_rates(k_Fe54_to_n_Fe53_derived)*Y(Fe54)) +
        (-screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(Fe54)*Y(N)*state.rho + screened_rates(k_Fe55_to_n_Fe54_derived)*Y(Fe55)) +
        (-screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(Fe54)*Y(H1)*state.rho + screened_rates(k_Co55_to_p_Fe54_derived)*Y(Co55)) +
        (-screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_Ni58_to_He4_Fe54_derived)*Y(Ni58)) +
        (screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(H1)*state.rho) +
        (screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*Y(N)*state.rho + -screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*Y(He4)*state.rho);

    ydot_nuc(Fe55) =
        (screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(Fe54)*Y(N)*state.rho + -screened_rates(k_Fe55_to_n_Fe54_derived)*Y(Fe55)) +
        (-screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(Fe55)*Y(N)*state.rho + screened_rates(k_Fe56_to_n_Fe55_derived)*Y(Fe56)) +
        (-screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(Fe55)*Y(H1)*state.rho + screened_rates(k_Co56_to_p_Fe55_derived)*Y(Co56)) +
        (screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*Y(Mn52)*state.rho + -screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*Y(N)*state.rho + -screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*Y(He4)*state.rho) +
        (screened_rates(k_Co55_to_Fe55_weaktab)*Y(Co55) + -screened_rates(k_Fe55_to_Co55_weaktab)*Y(Fe55));

    ydot_nuc(Fe56) =
        (screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(Fe55)*Y(N)*state.rho + -screened_rates(k_Fe56_to_n_Fe55_derived)*Y(Fe56)) +
        (-screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(Fe56)*Y(H1)*state.rho + screened_rates(k_Co57_to_p_Fe56_derived)*Y(Co57)) +
        (screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_Co56_to_Fe56_weaktab)*Y(Co56) + -screened_rates(k_Fe56_to_Co56_weaktab)*Y(Fe56));

    ydot_nuc(Co55) =
        (screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_Co55_to_He4_Mn51_derived)*Y(Co55)) +
        (-screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55_derived)*Y(Ni56)) +
        (screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(Fe54)*Y(H1)*state.rho + -screened_rates(k_Co55_to_p_Fe54_derived)*Y(Co55)) +
        (-screened_rates(k_n_Co55_to_Co56_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_Co56_to_n_Co55_derived)*Y(Co56)) +
        (-screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_Cu59_to_He4_Co55_derived)*Y(Cu59)) +
        (-screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*Y(Mn52)*state.rho) +
        (-screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (-screened_rates(k_Co55_to_Fe55_weaktab)*Y(Co55) + screened_rates(k_Fe55_to_Co55_weaktab)*Y(Fe55));

    ydot_nuc(Co56) =
        (screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(He4)*Y(Mn52)*state.rho + -screened_rates(k_Co56_to_He4_Mn52_derived)*Y(Co56)) +
        (screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(Fe55)*Y(H1)*state.rho + -screened_rates(k_Co56_to_p_Fe55_derived)*Y(Co56)) +
        (screened_rates(k_n_Co55_to_Co56_reaclib)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_Co56_to_n_Co55_derived)*Y(Co56)) +
        (-screened_rates(k_n_Co56_to_Co57_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Co57_to_n_Co56_derived)*Y(Co57)) +
        (-screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(Co56)*Y(H1)*state.rho + screened_rates(k_Ni57_to_p_Co56_derived)*Y(Ni57)) +
        (screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*Y(He4)*state.rho + -screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + -screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_Co56_to_Fe56_weaktab)*Y(Co56) + screened_rates(k_Fe56_to_Co56_weaktab)*Y(Fe56)) +
        (screened_rates(k_Ni56_to_Co56_weaktab)*Y(Ni56) + -screened_rates(k_Co56_to_Ni56_weaktab)*Y(Co56));

    ydot_nuc(Co57) =
        (screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(Fe56)*Y(H1)*state.rho + -screened_rates(k_Co57_to_p_Fe56_derived)*Y(Co57)) +
        (screened_rates(k_n_Co56_to_Co57_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_Co57_to_n_Co56_derived)*Y(Co57)) +
        (-screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_Ni58_to_p_Co57_derived)*Y(Ni58)) +
        (-screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*Y(N)*state.rho + -screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*Y(H1)*state.rho) +
        (screened_rates(k_Ni57_to_Co57_weaktab)*Y(Ni57) + -screened_rates(k_Co57_to_Ni57_weaktab)*Y(Co57));

    ydot_nuc(Ni56) =
        (screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(Co55)*Y(H1)*state.rho + -screened_rates(k_Ni56_to_p_Co55_derived)*Y(Ni56)) +
        (-screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_Ni57_to_n_Ni56_derived)*Y(Ni57)) +
        (-screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(He4)*Y(Ni56)*state.rho + screened_rates(k_Zn60_to_He4_Ni56_derived)*Y(Zn60)) +
        (-screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*Y(He4)*state.rho) +
        (screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*Y(Ni56)*state.rho) +
        (-screened_rates(k_Ni56_to_Co56_weaktab)*Y(Ni56) + screened_rates(k_Co56_to_Ni56_weaktab)*Y(Co56));

    ydot_nuc(Ni57) =
        (screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(Fe53)*Y(He4)*state.rho + -screened_rates(k_Ni57_to_He4_Fe53_derived)*Y(Ni57)) +
        (screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(Co56)*Y(H1)*state.rho + -screened_rates(k_Ni57_to_p_Co56_derived)*Y(Ni57)) +
        (screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_Ni57_to_n_Ni56_derived)*Y(Ni57)) +
        (-screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(Ni57)*Y(N)*state.rho + screened_rates(k_Ni58_to_n_Ni57_derived)*Y(Ni58)) +
        (-screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*Y(N)*state.rho + screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*Y(N)*state.rho + screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*Y(N)*state.rho + -screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*Y(Ni57)*state.rho) +
        (-screened_rates(k_Ni57_to_Co57_weaktab)*Y(Ni57) + screened_rates(k_Co57_to_Ni57_weaktab)*Y(Co57));

    ydot_nuc(Ni58) =
        (screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(Fe54)*Y(He4)*state.rho + -screened_rates(k_Ni58_to_He4_Fe54_derived)*Y(Ni58)) +
        (screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_Ni58_to_p_Co57_derived)*Y(Ni58)) +
        (screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(Ni57)*Y(N)*state.rho + -screened_rates(k_Ni58_to_n_Ni57_derived)*Y(Ni58)) +
        (-screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(Ni58)*Y(H1)*state.rho + screened_rates(k_Cu59_to_p_Ni58_derived)*Y(Cu59)) +
        (screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (-screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*Y(N)*state.rho + screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*Y(He4)*state.rho);

    ydot_nuc(Cu59) =
        (screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_Cu59_to_He4_Co55_derived)*Y(Cu59)) +
        (screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(Ni58)*Y(H1)*state.rho + -screened_rates(k_Cu59_to_p_Ni58_derived)*Y(Cu59)) +
        (-screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_Zn60_to_p_Cu59_derived)*Y(Zn60)) +
        (-screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*Y(Ni56)*state.rho);

    ydot_nuc(Zn60) =
        (screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(He4)*Y(Ni56)*state.rho + -screened_rates(k_Zn60_to_He4_Ni56_derived)*Y(Zn60)) +
        (screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_Zn60_to_p_Cu59_derived)*Y(Zn60)) +
        (-screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*Y(N)*state.rho + screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*Y(Ni57)*state.rho);

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
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates)
{

    amrex::Real scratch;

    scratch = -screened_rates(k_n_Co55_to_Co56_reaclib)*Y(Co55)*state.rho - screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*state.rho - screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*state.rho - screened_rates(k_n_Co56_to_Co57_reaclib)*Y(Co56)*state.rho - screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*state.rho - screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*state.rho - screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(Fe52)*state.rho - screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*state.rho - screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(Fe53)*state.rho - screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(Fe54)*state.rho - screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(Fe55)*state.rho - screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(Mn50)*state.rho - screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(Mn51)*state.rho - screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*state.rho - screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(Ni56)*state.rho - screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*state.rho - screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*state.rho - screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(Ni57)*state.rho - screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*state.rho - screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*state.rho - screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*state.rho - screened_rates(k_n_to_p_weaktab);
    jac.set(N, N, scratch);

    scratch = screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho + screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*state.rho + screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*state.rho + screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho + screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*state.rho + screened_rates(k_p_to_n_weaktab);
    jac.set(N, H1, scratch);

    scratch = screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*state.rho + screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*state.rho + screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*state.rho + screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*state.rho + screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(Mn52)*state.rho + screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(Ni57)*state.rho;
    jac.set(N, He4, scratch);

    jac.set(N, C12, 0.0);

    jac.set(N, C13, 0.0);

    jac.set(N, C14, 0.0);

    jac.set(N, N13, 0.0);

    jac.set(N, N14, 0.0);

    jac.set(N, N15, 0.0);

    jac.set(N, O16, 0.0);

    jac.set(N, O17, 0.0);

    jac.set(N, O18, 0.0);

    jac.set(N, F18, 0.0);

    jac.set(N, F19, 0.0);

    jac.set(N, Ne20, 0.0);

    jac.set(N, Ne21, 0.0);

    jac.set(N, Ne22, 0.0);

    jac.set(N, Na22, 0.0);

    jac.set(N, Na23, 0.0);

    jac.set(N, Mg24, 0.0);

    jac.set(N, Mg25, 0.0);

    jac.set(N, Al26, 0.0);

    jac.set(N, Al27, 0.0);

    jac.set(N, Si28, 0.0);

    jac.set(N, Si29, 0.0);

    jac.set(N, P30, 0.0);

    jac.set(N, P31, 0.0);

    jac.set(N, S32, 0.0);

    jac.set(N, S33, 0.0);

    jac.set(N, Cl34, 0.0);

    jac.set(N, Cl35, 0.0);

    jac.set(N, Ar36, 0.0);

    jac.set(N, Ar37, 0.0);

    jac.set(N, K38, 0.0);

    jac.set(N, K39, 0.0);

    jac.set(N, Ca40, 0.0);

    jac.set(N, Ca41, 0.0);

    jac.set(N, Sc42, 0.0);

    jac.set(N, Sc43, 0.0);

    jac.set(N, Ti44, 0.0);

    jac.set(N, Ti45, 0.0);

    jac.set(N, V46, 0.0);

    jac.set(N, V47, 0.0);

    jac.set(N, Cr48, 0.0);

    jac.set(N, Cr49, 0.0);

    scratch = -screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(N)*state.rho;
    jac.set(N, Mn50, scratch);

    scratch = screened_rates(k_Mn51_to_n_Mn50_derived) - screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(N)*state.rho;
    jac.set(N, Mn51, scratch);

    scratch = screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*state.rho + screened_rates(k_Mn52_to_n_Mn51_derived) + screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(H1)*state.rho;
    jac.set(N, Mn52, scratch);

    scratch = -screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(N)*state.rho - screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(N)*state.rho;
    jac.set(N, Fe52, scratch);

    scratch = screened_rates(k_Fe53_to_n_Fe52_derived) + screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(He4)*state.rho - screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(N)*state.rho;
    jac.set(N, Fe53, scratch);

    scratch = screened_rates(k_Fe54_to_n_Fe53_derived) + screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(He4)*state.rho - screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(N)*state.rho;
    jac.set(N, Fe54, scratch);

    scratch = screened_rates(k_Fe55_to_n_Fe54_derived) + screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(He4)*state.rho - screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(N)*state.rho + screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(H1)*state.rho;
    jac.set(N, Fe55, scratch);

    scratch = screened_rates(k_Fe56_to_n_Fe55_derived) + screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(H1)*state.rho;
    jac.set(N, Fe56, scratch);

    scratch = -screened_rates(k_n_Co55_to_Co56_reaclib)*Y(N)*state.rho - screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(N)*state.rho - screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(N)*state.rho;
    jac.set(N, Co55, scratch);

    scratch = screened_rates(k_Co56_to_n_Co55_derived) + screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(He4)*state.rho - screened_rates(k_n_Co56_to_Co57_reaclib)*Y(N)*state.rho - screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(N)*state.rho + screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(H1)*state.rho;
    jac.set(N, Co56, scratch);

    scratch = screened_rates(k_Co57_to_n_Co56_derived) + screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(H1)*state.rho;
    jac.set(N, Co57, scratch);

    scratch = -screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(N)*state.rho;
    jac.set(N, Ni56, scratch);

    scratch = screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*state.rho + screened_rates(k_Ni57_to_n_Ni56_derived) - screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(N)*state.rho;
    jac.set(N, Ni57, scratch);

    scratch = screened_rates(k_Ni58_to_n_Ni57_derived) - screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(N)*state.rho;
    jac.set(N, Ni58, scratch);

    scratch = -screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(N)*state.rho;
    jac.set(N, Cu59, scratch);

    scratch = -screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(N)*state.rho;
    jac.set(N, Zn60, scratch);

    scratch = screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*state.rho + screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*state.rho + screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*state.rho + screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*state.rho + screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*state.rho + screened_rates(k_n_to_p_weaktab);
    jac.set(H1, N, scratch);

    scratch = -screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*state.rho - screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*state.rho - screened_rates(k_p_Ar37_to_K38_reaclib)*Y(Ar37)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho - screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*state.rho - screened_rates(k_p_C14_to_N15_reaclib)*Y(C14)*state.rho - screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*state.rho - screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(Ca41)*state.rho - screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(Cl35)*state.rho - screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*state.rho - screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(Co55)*state.rho - screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*state.rho - screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(Co56)*state.rho - screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho - screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*state.rho - screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(Co57)*state.rho - screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*state.rho - screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*state.rho - screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(Cr49)*state.rho - screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho - screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(Cu59)*state.rho - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho - screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*state.rho - screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(Fe54)*state.rho - screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho - screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(Fe55)*state.rho - screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*state.rho - screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*state.rho - screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(Fe56)*state.rho - screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho - screened_rates(k_p_K39_to_Ca40_reaclib)*Y(K39)*state.rho - screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*state.rho - screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(Mg25)*state.rho - screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*state.rho - screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(Mn51)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho - screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(Mn52)*state.rho - screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*state.rho - screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*state.rho - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho - screened_rates(k_p_Na23_to_C12_C12_derived)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*state.rho - screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*state.rho - screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*state.rho - screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(Ne22)*state.rho - screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(Ni58)*state.rho - screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho - screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho - screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*state.rho - screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho - screened_rates(k_p_P31_to_O16_O16_derived)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*state.rho - screened_rates(k_p_S33_to_Cl34_reaclib)*Y(S33)*state.rho - screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*state.rho - screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho - screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(Sc43)*state.rho - screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*state.rho - screened_rates(k_p_Si29_to_P30_reaclib)*Y(Si29)*state.rho - screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*state.rho - screened_rates(k_p_Ti45_to_V46_reaclib)*Y(Ti45)*state.rho - screened_rates(k_p_V47_to_Cr48_reaclib)*Y(V47)*state.rho - screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*state.rho - screened_rates(k_p_to_n_weaktab);
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*state.rho + screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho + screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*state.rho + screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*state.rho + screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*state.rho + screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*state.rho + screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*state.rho + screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho + screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*state.rho + screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*state.rho + screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*state.rho + screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*state.rho + screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(K38)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(Mn50)*state.rho + screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(Mn51)*state.rho + screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(Mn52)*state.rho + screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho + screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(Na22)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho + screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(Ni56)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho + screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(P30)*state.rho + screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(S32)*state.rho + screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(Sc42)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho + screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(Ti44)*state.rho + screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(V46)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho + screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(Ne20)*state.rho + screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_p_C13_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, C13, scratch);

    scratch = -screened_rates(k_p_C14_to_N15_reaclib)*Y(H1)*state.rho;
    jac.set(H1, C14, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho + screened_rates(k_N13_to_p_C12_derived);
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho + screened_rates(k_N14_to_p_C13_derived);
    jac.set(H1, N14, scratch);

    scratch = screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho + screened_rates(k_N15_to_p_C14_derived) - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(H1, N15, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho + screened_rates(k_O16_to_p_N15_derived) - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = -screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O17, scratch);

    scratch = -screened_rates(k_p_O18_to_F19_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O18, scratch);

    scratch = screened_rates(k_F18_to_p_O17_derived) + screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18_derived) + screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(He4)*state.rho - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(H1, F19, scratch);

    scratch = screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho + screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(H1, Ne20, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ne21, scratch);

    scratch = -screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ne22, scratch);

    scratch = screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*state.rho + screened_rates(k_Na22_to_p_Ne21_derived);
    jac.set(H1, Na22, scratch);

    scratch = screened_rates(k_Na23_to_p_Ne22_derived) - screened_rates(k_p_Na23_to_C12_C12_derived)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Na23, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(H1, Mg24, scratch);

    scratch = -screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(H1)*state.rho;
    jac.set(H1, Mg25, scratch);

    scratch = screened_rates(k_Al26_to_p_Mg25_derived) + screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(He4)*state.rho;
    jac.set(H1, Al26, scratch);

    scratch = -screened_rates(k_p_Al27_to_C12_O16_derived)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho + screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(H1)*state.rho - screened_rates(k_p_Si29_to_P30_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Si29, scratch);

    scratch = screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*state.rho + screened_rates(k_P30_to_p_Si29_derived);
    jac.set(H1, P30, scratch);

    scratch = -screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(H1)*state.rho - screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho - screened_rates(k_p_P31_to_O16_O16_derived)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(H1, P31, scratch);

    scratch = screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*state.rho + screened_rates(k_S32_to_p_P31_derived);
    jac.set(H1, S32, scratch);

    scratch = -screened_rates(k_p_S33_to_Cl34_reaclib)*Y(H1)*state.rho - screened_rates(k_p_S33_to_He4_P30_derived)*Y(H1)*state.rho;
    jac.set(H1, S33, scratch);

    scratch = screened_rates(k_Cl34_to_p_S33_derived) + screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(He4)*state.rho;
    jac.set(H1, Cl34, scratch);

    scratch = -screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_p_Cl35_derived) + screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(He4)*state.rho;
    jac.set(H1, Ar36, scratch);

    scratch = -screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(H1)*state.rho - screened_rates(k_p_Ar37_to_K38_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ar37, scratch);

    scratch = screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*state.rho + screened_rates(k_K38_to_p_Ar37_derived);
    jac.set(H1, K38, scratch);

    scratch = -screened_rates(k_p_K39_to_Ca40_reaclib)*Y(H1)*state.rho - screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(H1)*state.rho;
    jac.set(H1, K39, scratch);

    scratch = screened_rates(k_Ca40_to_p_K39_derived) + screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(He4)*state.rho;
    jac.set(H1, Ca40, scratch);

    scratch = -screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(H1)*state.rho - screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ca41, scratch);

    scratch = screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*state.rho + screened_rates(k_Sc42_to_p_Ca41_derived);
    jac.set(H1, Sc42, scratch);

    scratch = -screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Sc43, scratch);

    scratch = screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*state.rho + screened_rates(k_Ti44_to_p_Sc43_derived);
    jac.set(H1, Ti44, scratch);

    scratch = -screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(H1)*state.rho - screened_rates(k_p_Ti45_to_V46_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ti45, scratch);

    scratch = screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*state.rho + screened_rates(k_V46_to_p_Ti45_derived);
    jac.set(H1, V46, scratch);

    scratch = -screened_rates(k_p_V47_to_Cr48_reaclib)*Y(H1)*state.rho - screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(H1)*state.rho;
    jac.set(H1, V47, scratch);

    scratch = screened_rates(k_Cr48_to_p_V47_derived) + screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(He4)*state.rho;
    jac.set(H1, Cr48, scratch);

    scratch = screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(H1)*state.rho - screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Cr49, scratch);

    scratch = screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*state.rho + screened_rates(k_Mn50_to_p_Cr49_derived);
    jac.set(H1, Mn50, scratch);

    scratch = screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(H1)*state.rho;
    jac.set(H1, Mn51, scratch);

    scratch = screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(H1)*state.rho - screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(H1)*state.rho;
    jac.set(H1, Mn52, scratch);

    scratch = screened_rates(k_Fe52_to_p_Mn51_derived) + screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(He4)*state.rho + screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(N)*state.rho;
    jac.set(H1, Fe52, scratch);

    scratch = screened_rates(k_Fe53_to_p_Mn52_derived) + screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(H1)*state.rho;
    jac.set(H1, Fe53, scratch);

    scratch = screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(He4)*state.rho - screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(H1)*state.rho;
    jac.set(H1, Fe54, scratch);

    scratch = -screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(H1)*state.rho - screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(H1)*state.rho;
    jac.set(H1, Fe55, scratch);

    scratch = -screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(H1)*state.rho;
    jac.set(H1, Fe56, scratch);

    scratch = screened_rates(k_Co55_to_p_Fe54_derived) + screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(He4)*state.rho + screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(N)*state.rho - screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho - screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Co55, scratch);

    scratch = screened_rates(k_Co56_to_p_Fe55_derived) + screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(N)*state.rho - screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(H1)*state.rho - screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(H1)*state.rho;
    jac.set(H1, Co56, scratch);

    scratch = screened_rates(k_Co57_to_p_Fe56_derived) - screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(H1)*state.rho;
    jac.set(H1, Co57, scratch);

    scratch = screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*state.rho + screened_rates(k_Ni56_to_p_Co55_derived) + screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(N)*state.rho;
    jac.set(H1, Ni56, scratch);

    scratch = screened_rates(k_Ni57_to_p_Co56_derived) + screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(N)*state.rho;
    jac.set(H1, Ni57, scratch);

    scratch = screened_rates(k_Ni58_to_p_Co57_derived) - screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(H1)*state.rho;
    jac.set(H1, Ni58, scratch);

    scratch = screened_rates(k_Cu59_to_p_Ni58_derived) - screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_p_Cu59_derived);
    jac.set(H1, Zn60, scratch);

    scratch = screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*state.rho + screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*state.rho + screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*state.rho + screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*state.rho + screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*state.rho + screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*state.rho;
    jac.set(He4, N, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho + screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*state.rho + screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*state.rho + screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho + screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*state.rho + screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*state.rho + screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*state.rho + screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho + screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*state.rho + screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho + screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*state.rho + screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*state.rho + screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*state.rho + screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho + screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho + screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho + screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*state.rho + screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho + screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho + screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho + screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho + screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*state.rho + screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho + screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*state.rho + screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*state.rho + screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_He4_Al26_to_P30_reaclib)*Y(Al26)*state.rho - screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*state.rho - screened_rates(k_He4_Al27_to_P31_reaclib)*Y(Al27)*state.rho - screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(Ar36)*state.rho - screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*state.rho - screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(Ar37)*state.rho - screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho - screened_rates(k_He4_C14_to_O18_reaclib)*Y(C14)*state.rho - screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(Ca40)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*state.rho - screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(Ca41)*state.rho - screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(Cl34)*state.rho - screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*state.rho - screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(Cl35)*state.rho - screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(Co55)*state.rho - screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*state.rho - screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*state.rho - screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(Cr48)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*state.rho - screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*state.rho - screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*state.rho - screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*state.rho - screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(Fe53)*state.rho - screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*state.rho - screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*state.rho - screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(Fe54)*state.rho - screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*state.rho - screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*state.rho - screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) - screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(K38)*state.rho - screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(K38)*state.rho - screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(K39)*state.rho - screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(Mg25)*state.rho - screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(Mn50)*state.rho - screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(Mn51)*state.rho - screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(Mn51)*state.rho - screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(Mn52)*state.rho - screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(Mn52)*state.rho - screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(Mn52)*state.rho - screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho - screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho - screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(Na22)*state.rho - screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(Na22)*state.rho - screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(Na23)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(Ne21)*state.rho - screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(Ni56)*state.rho - screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(Ni56)*state.rho - screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(Ni57)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho - screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(O17)*state.rho - screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(O18)*state.rho - screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(P30)*state.rho - screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(P30)*state.rho - screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(P31)*state.rho - screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(S32)*state.rho - screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(S32)*state.rho - screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(S33)*state.rho - screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(Sc42)*state.rho - screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(Sc42)*state.rho - screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(Sc43)*state.rho - screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_S32_reaclib)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho - screened_rates(k_He4_Si29_to_S33_reaclib)*Y(Si29)*state.rho - screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(Ti44)*state.rho - screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(Ti44)*state.rho - screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(Ti45)*state.rho - screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(V46)*state.rho - screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(V46)*state.rho - screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(V47)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho + screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(Ne20)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho + 3.0*screened_rates(k_C12_to_He4_He4_He4_derived) - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    jac.set(He4, C13, 0.0);

    scratch = -screened_rates(k_He4_C14_to_O18_reaclib)*Y(He4)*state.rho;
    jac.set(He4, C14, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = -screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho;
    jac.set(He4, N15, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho + screened_rates(k_O16_to_He4_C12_derived) + screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = -screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*state.rho + screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O17, scratch);

    scratch = -screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(He4)*state.rho + screened_rates(k_O18_to_He4_C14_derived) + screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O18, scratch);

    scratch = screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_He4_F18_to_Na22_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = screened_rates(k_F19_to_He4_N15_derived) - screened_rates(k_He4_F19_to_Na23_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(He4)*state.rho + screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho;
    jac.set(He4, F19, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(He4, Ne20, scratch);

    scratch = -screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne21_to_He4_O17_derived) + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho;
    jac.set(He4, Ne21, scratch);

    scratch = screened_rates(k_Ne22_to_He4_O18_derived) + screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(H1)*state.rho;
    jac.set(He4, Ne22, scratch);

    scratch = -screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18_derived);
    jac.set(He4, Na22, scratch);

    scratch = -screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(He4)*state.rho + screened_rates(k_Na23_to_He4_F19_derived) + screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Na23, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(He4, Mg24, scratch);

    scratch = -screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(He4)*state.rho + screened_rates(k_Mg25_to_He4_Ne21_derived) + screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(H1)*state.rho;
    jac.set(He4, Mg25, scratch);

    scratch = screened_rates(k_Al26_to_He4_Na22_derived) - screened_rates(k_He4_Al26_to_P30_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Al26, scratch);

    scratch = screened_rates(k_Al27_to_He4_Na23_derived) - screened_rates(k_He4_Al27_to_P31_reaclib)*Y(He4)*state.rho + screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived);
    jac.set(He4, Si28, scratch);

    scratch = -screened_rates(k_He4_Si29_to_S33_reaclib)*Y(He4)*state.rho + screened_rates(k_Si29_to_He4_Mg25_derived) + screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(H1)*state.rho;
    jac.set(He4, Si29, scratch);

    scratch = -screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*state.rho + screened_rates(k_P30_to_He4_Al26_derived);
    jac.set(He4, P30, scratch);

    scratch = -screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27_derived) + screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(He4, P31, scratch);

    scratch = -screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*state.rho + screened_rates(k_S32_to_He4_Si28_derived);
    jac.set(He4, S32, scratch);

    scratch = -screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(He4)*state.rho + screened_rates(k_S33_to_He4_Si29_derived) + screened_rates(k_p_S33_to_He4_P30_derived)*Y(H1)*state.rho;
    jac.set(He4, S33, scratch);

    scratch = screened_rates(k_Cl34_to_He4_P30_derived) - screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Cl34, scratch);

    scratch = screened_rates(k_Cl35_to_He4_P31_derived) - screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(He4)*state.rho + screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_He4_S32_derived) - screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(He4)*state.rho;
    jac.set(He4, Ar36, scratch);

    scratch = screened_rates(k_Ar37_to_He4_S33_derived) - screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(He4)*state.rho + screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(H1)*state.rho;
    jac.set(He4, Ar37, scratch);

    scratch = -screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*state.rho + screened_rates(k_K38_to_He4_Cl34_derived);
    jac.set(He4, K38, scratch);

    scratch = -screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35_derived) + screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(H1)*state.rho;
    jac.set(He4, K39, scratch);

    scratch = screened_rates(k_Ca40_to_He4_Ar36_derived) - screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(He4)*state.rho;
    jac.set(He4, Ca40, scratch);

    scratch = screened_rates(k_Ca41_to_He4_Ar37_derived) - screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(He4)*state.rho + screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(H1)*state.rho;
    jac.set(He4, Ca41, scratch);

    scratch = -screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*state.rho + screened_rates(k_Sc42_to_He4_K38_derived);
    jac.set(He4, Sc42, scratch);

    scratch = -screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(He4)*state.rho + screened_rates(k_Sc43_to_He4_K39_derived) + screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Sc43, scratch);

    scratch = -screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40_derived);
    jac.set(He4, Ti44, scratch);

    scratch = -screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(He4)*state.rho + screened_rates(k_Ti45_to_He4_Ca41_derived) + screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(H1)*state.rho;
    jac.set(He4, Ti45, scratch);

    scratch = -screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*state.rho + screened_rates(k_V46_to_He4_Sc42_derived);
    jac.set(He4, V46, scratch);

    scratch = -screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(He4)*state.rho + screened_rates(k_V47_to_He4_Sc43_derived) + screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(H1)*state.rho;
    jac.set(He4, V47, scratch);

    scratch = screened_rates(k_Cr48_to_He4_Ti44_derived) - screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Cr48, scratch);

    scratch = screened_rates(k_Cr49_to_He4_Ti45_derived) - screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(He4)*state.rho + screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(H1)*state.rho;
    jac.set(He4, Cr49, scratch);

    scratch = -screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*state.rho + screened_rates(k_Mn50_to_He4_V46_derived);
    jac.set(He4, Mn50, scratch);

    scratch = -screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*state.rho + screened_rates(k_Mn51_to_He4_V47_derived) + screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(H1)*state.rho;
    jac.set(He4, Mn51, scratch);

    scratch = -screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*state.rho - screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*state.rho + screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(H1)*state.rho;
    jac.set(He4, Mn52, scratch);

    scratch = screened_rates(k_Fe52_to_He4_Cr48_derived) - screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Fe52, scratch);

    scratch = -screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(He4)*state.rho - screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(He4)*state.rho + screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(H1)*state.rho;
    jac.set(He4, Fe53, scratch);

    scratch = -screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(He4)*state.rho - screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(He4)*state.rho + screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(H1)*state.rho;
    jac.set(He4, Fe54, scratch);

    scratch = -screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(He4)*state.rho + screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(H1)*state.rho;
    jac.set(He4, Fe55, scratch);

    jac.set(He4, Fe56, 0.0);

    scratch = screened_rates(k_Co55_to_He4_Mn51_derived) - screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(He4)*state.rho + screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(N)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho;
    jac.set(He4, Co55, scratch);

    scratch = screened_rates(k_Co56_to_He4_Mn52_derived) - screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(He4)*state.rho + screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(H1)*state.rho;
    jac.set(He4, Co56, scratch);

    scratch = screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Co57, scratch);

    scratch = -screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_derived) + screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(N)*state.rho;
    jac.set(He4, Ni56, scratch);

    scratch = -screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*state.rho + screened_rates(k_Ni57_to_He4_Fe53_derived) + screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(N)*state.rho;
    jac.set(He4, Ni57, scratch);

    scratch = screened_rates(k_Ni58_to_He4_Fe54_derived) + screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(N)*state.rho + screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(H1)*state.rho;
    jac.set(He4, Ni58, scratch);

    scratch = screened_rates(k_Cu59_to_He4_Co55_derived) + screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(N)*state.rho + screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_He4_Ni56_derived) + screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(N)*state.rho;
    jac.set(He4, Zn60, scratch);

    jac.set(C12, N, 0.0);

    scratch = screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho + 2.0*screened_rates(k_p_Na23_to_C12_C12_derived)*Y(Na23)*state.rho + screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(Si28)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(Ne20)*state.rho - screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(Ne20)*state.rho - screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho - screened_rates(k_C12_to_He4_He4_He4_derived) - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    jac.set(C12, C13, 0.0);

    jac.set(C12, C14, 0.0);

    scratch = screened_rates(k_N13_to_p_C12_derived);
    jac.set(C12, N13, scratch);

    jac.set(C12, N14, 0.0);

    scratch = screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho;
    jac.set(C12, N15, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho + screened_rates(k_O16_to_He4_C12_derived);
    jac.set(C12, O16, scratch);

    jac.set(C12, O17, 0.0);

    jac.set(C12, O18, 0.0);

    jac.set(C12, F18, 0.0);

    jac.set(C12, F19, 0.0);

    scratch = -screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*state.rho;
    jac.set(C12, Ne20, scratch);

    jac.set(C12, Ne21, 0.0);

    jac.set(C12, Ne22, 0.0);

    jac.set(C12, Na22, 0.0);

    scratch = 2.0*screened_rates(k_p_Na23_to_C12_C12_derived)*Y(H1)*state.rho;
    jac.set(C12, Na23, scratch);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho;
    jac.set(C12, Mg24, scratch);

    jac.set(C12, Mg25, 0.0);

    jac.set(C12, Al26, 0.0);

    scratch = screened_rates(k_p_Al27_to_C12_O16_derived)*Y(H1)*state.rho;
    jac.set(C12, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*state.rho;
    jac.set(C12, Si28, scratch);

    jac.set(C12, Si29, 0.0);

    jac.set(C12, P30, 0.0);

    scratch = screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(H1)*state.rho;
    jac.set(C12, P31, scratch);

    jac.set(C12, S32, 0.0);

    jac.set(C12, S33, 0.0);

    jac.set(C12, Cl34, 0.0);

    jac.set(C12, Cl35, 0.0);

    jac.set(C12, Ar36, 0.0);

    jac.set(C12, Ar37, 0.0);

    jac.set(C12, K38, 0.0);

    jac.set(C12, K39, 0.0);

    jac.set(C12, Ca40, 0.0);

    jac.set(C12, Ca41, 0.0);

    jac.set(C12, Sc42, 0.0);

    jac.set(C12, Sc43, 0.0);

    jac.set(C12, Ti44, 0.0);

    jac.set(C12, Ti45, 0.0);

    jac.set(C12, V46, 0.0);

    jac.set(C12, V47, 0.0);

    jac.set(C12, Cr48, 0.0);

    jac.set(C12, Cr49, 0.0);

    jac.set(C12, Mn50, 0.0);

    jac.set(C12, Mn51, 0.0);

    jac.set(C12, Mn52, 0.0);

    jac.set(C12, Fe52, 0.0);

    jac.set(C12, Fe53, 0.0);

    jac.set(C12, Fe54, 0.0);

    jac.set(C12, Fe55, 0.0);

    jac.set(C12, Fe56, 0.0);

    jac.set(C12, Co55, 0.0);

    jac.set(C12, Co56, 0.0);

    jac.set(C12, Co57, 0.0);

    jac.set(C12, Ni56, 0.0);

    jac.set(C12, Ni57, 0.0);

    jac.set(C12, Ni58, 0.0);

    jac.set(C12, Cu59, 0.0);

    jac.set(C12, Zn60, 0.0);

    jac.set(C13, N, 0.0);

    scratch = -screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*state.rho;
    jac.set(C13, H1, scratch);

    jac.set(C13, He4, 0.0);

    jac.set(C13, C12, 0.0);

    scratch = -screened_rates(k_p_C13_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(C13, C13, scratch);

    jac.set(C13, C14, 0.0);

    scratch = screened_rates(k_N13_to_C13_reaclib);
    jac.set(C13, N13, scratch);

    scratch = screened_rates(k_N14_to_p_C13_derived);
    jac.set(C13, N14, scratch);

    jac.set(C13, N15, 0.0);

    jac.set(C13, O16, 0.0);

    jac.set(C13, O17, 0.0);

    jac.set(C13, O18, 0.0);

    jac.set(C13, F18, 0.0);

    jac.set(C13, F19, 0.0);

    jac.set(C13, Ne20, 0.0);

    jac.set(C13, Ne21, 0.0);

    jac.set(C13, Ne22, 0.0);

    jac.set(C13, Na22, 0.0);

    jac.set(C13, Na23, 0.0);

    jac.set(C13, Mg24, 0.0);

    jac.set(C13, Mg25, 0.0);

    jac.set(C13, Al26, 0.0);

    jac.set(C13, Al27, 0.0);

    jac.set(C13, Si28, 0.0);

    jac.set(C13, Si29, 0.0);

    jac.set(C13, P30, 0.0);

    jac.set(C13, P31, 0.0);

    jac.set(C13, S32, 0.0);

    jac.set(C13, S33, 0.0);

    jac.set(C13, Cl34, 0.0);

    jac.set(C13, Cl35, 0.0);

    jac.set(C13, Ar36, 0.0);

    jac.set(C13, Ar37, 0.0);

    jac.set(C13, K38, 0.0);

    jac.set(C13, K39, 0.0);

    jac.set(C13, Ca40, 0.0);

    jac.set(C13, Ca41, 0.0);

    jac.set(C13, Sc42, 0.0);

    jac.set(C13, Sc43, 0.0);

    jac.set(C13, Ti44, 0.0);

    jac.set(C13, Ti45, 0.0);

    jac.set(C13, V46, 0.0);

    jac.set(C13, V47, 0.0);

    jac.set(C13, Cr48, 0.0);

    jac.set(C13, Cr49, 0.0);

    jac.set(C13, Mn50, 0.0);

    jac.set(C13, Mn51, 0.0);

    jac.set(C13, Mn52, 0.0);

    jac.set(C13, Fe52, 0.0);

    jac.set(C13, Fe53, 0.0);

    jac.set(C13, Fe54, 0.0);

    jac.set(C13, Fe55, 0.0);

    jac.set(C13, Fe56, 0.0);

    jac.set(C13, Co55, 0.0);

    jac.set(C13, Co56, 0.0);

    jac.set(C13, Co57, 0.0);

    jac.set(C13, Ni56, 0.0);

    jac.set(C13, Ni57, 0.0);

    jac.set(C13, Ni58, 0.0);

    jac.set(C13, Cu59, 0.0);

    jac.set(C13, Zn60, 0.0);

    jac.set(C14, N, 0.0);

    scratch = -screened_rates(k_p_C14_to_N15_reaclib)*Y(C14)*state.rho;
    jac.set(C14, H1, scratch);

    scratch = -screened_rates(k_He4_C14_to_O18_reaclib)*Y(C14)*state.rho;
    jac.set(C14, He4, scratch);

    jac.set(C14, C12, 0.0);

    jac.set(C14, C13, 0.0);

    scratch = -screened_rates(k_C14_to_N14_reaclib) - screened_rates(k_He4_C14_to_O18_reaclib)*Y(He4)*state.rho - screened_rates(k_p_C14_to_N15_reaclib)*Y(H1)*state.rho;
    jac.set(C14, C14, scratch);

    jac.set(C14, N13, 0.0);

    jac.set(C14, N14, 0.0);

    scratch = screened_rates(k_N15_to_p_C14_derived);
    jac.set(C14, N15, scratch);

    jac.set(C14, O16, 0.0);

    jac.set(C14, O17, 0.0);

    scratch = screened_rates(k_O18_to_He4_C14_derived);
    jac.set(C14, O18, scratch);

    jac.set(C14, F18, 0.0);

    jac.set(C14, F19, 0.0);

    jac.set(C14, Ne20, 0.0);

    jac.set(C14, Ne21, 0.0);

    jac.set(C14, Ne22, 0.0);

    jac.set(C14, Na22, 0.0);

    jac.set(C14, Na23, 0.0);

    jac.set(C14, Mg24, 0.0);

    jac.set(C14, Mg25, 0.0);

    jac.set(C14, Al26, 0.0);

    jac.set(C14, Al27, 0.0);

    jac.set(C14, Si28, 0.0);

    jac.set(C14, Si29, 0.0);

    jac.set(C14, P30, 0.0);

    jac.set(C14, P31, 0.0);

    jac.set(C14, S32, 0.0);

    jac.set(C14, S33, 0.0);

    jac.set(C14, Cl34, 0.0);

    jac.set(C14, Cl35, 0.0);

    jac.set(C14, Ar36, 0.0);

    jac.set(C14, Ar37, 0.0);

    jac.set(C14, K38, 0.0);

    jac.set(C14, K39, 0.0);

    jac.set(C14, Ca40, 0.0);

    jac.set(C14, Ca41, 0.0);

    jac.set(C14, Sc42, 0.0);

    jac.set(C14, Sc43, 0.0);

    jac.set(C14, Ti44, 0.0);

    jac.set(C14, Ti45, 0.0);

    jac.set(C14, V46, 0.0);

    jac.set(C14, V47, 0.0);

    jac.set(C14, Cr48, 0.0);

    jac.set(C14, Cr49, 0.0);

    jac.set(C14, Mn50, 0.0);

    jac.set(C14, Mn51, 0.0);

    jac.set(C14, Mn52, 0.0);

    jac.set(C14, Fe52, 0.0);

    jac.set(C14, Fe53, 0.0);

    jac.set(C14, Fe54, 0.0);

    jac.set(C14, Fe55, 0.0);

    jac.set(C14, Fe56, 0.0);

    jac.set(C14, Co55, 0.0);

    jac.set(C14, Co56, 0.0);

    jac.set(C14, Co57, 0.0);

    jac.set(C14, Ni56, 0.0);

    jac.set(C14, Ni57, 0.0);

    jac.set(C14, Ni58, 0.0);

    jac.set(C14, Cu59, 0.0);

    jac.set(C14, Zn60, 0.0);

    jac.set(N13, N, 0.0);

    scratch = screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    jac.set(N13, C13, 0.0);

    jac.set(N13, C14, 0.0);

    scratch = -screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_N13_to_C13_reaclib) - screened_rates(k_N13_to_p_C12_derived);
    jac.set(N13, N13, scratch);

    jac.set(N13, N14, 0.0);

    jac.set(N13, N15, 0.0);

    scratch = screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(N13, O16, scratch);

    jac.set(N13, O17, 0.0);

    jac.set(N13, O18, 0.0);

    jac.set(N13, F18, 0.0);

    jac.set(N13, F19, 0.0);

    jac.set(N13, Ne20, 0.0);

    jac.set(N13, Ne21, 0.0);

    jac.set(N13, Ne22, 0.0);

    jac.set(N13, Na22, 0.0);

    jac.set(N13, Na23, 0.0);

    jac.set(N13, Mg24, 0.0);

    jac.set(N13, Mg25, 0.0);

    jac.set(N13, Al26, 0.0);

    jac.set(N13, Al27, 0.0);

    jac.set(N13, Si28, 0.0);

    jac.set(N13, Si29, 0.0);

    jac.set(N13, P30, 0.0);

    jac.set(N13, P31, 0.0);

    jac.set(N13, S32, 0.0);

    jac.set(N13, S33, 0.0);

    jac.set(N13, Cl34, 0.0);

    jac.set(N13, Cl35, 0.0);

    jac.set(N13, Ar36, 0.0);

    jac.set(N13, Ar37, 0.0);

    jac.set(N13, K38, 0.0);

    jac.set(N13, K39, 0.0);

    jac.set(N13, Ca40, 0.0);

    jac.set(N13, Ca41, 0.0);

    jac.set(N13, Sc42, 0.0);

    jac.set(N13, Sc43, 0.0);

    jac.set(N13, Ti44, 0.0);

    jac.set(N13, Ti45, 0.0);

    jac.set(N13, V46, 0.0);

    jac.set(N13, V47, 0.0);

    jac.set(N13, Cr48, 0.0);

    jac.set(N13, Cr49, 0.0);

    jac.set(N13, Mn50, 0.0);

    jac.set(N13, Mn51, 0.0);

    jac.set(N13, Mn52, 0.0);

    jac.set(N13, Fe52, 0.0);

    jac.set(N13, Fe53, 0.0);

    jac.set(N13, Fe54, 0.0);

    jac.set(N13, Fe55, 0.0);

    jac.set(N13, Fe56, 0.0);

    jac.set(N13, Co55, 0.0);

    jac.set(N13, Co56, 0.0);

    jac.set(N13, Co57, 0.0);

    jac.set(N13, Ni56, 0.0);

    jac.set(N13, Ni57, 0.0);

    jac.set(N13, Ni58, 0.0);

    jac.set(N13, Cu59, 0.0);

    jac.set(N13, Zn60, 0.0);

    jac.set(N14, N, 0.0);

    scratch = screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*state.rho + screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho;
    jac.set(N14, H1, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    jac.set(N14, C12, 0.0);

    scratch = screened_rates(k_p_C13_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(N14, C13, scratch);

    scratch = screened_rates(k_C14_to_N14_reaclib);
    jac.set(N14, C14, scratch);

    jac.set(N14, N13, 0.0);

    scratch = -screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho - screened_rates(k_N14_to_p_C13_derived);
    jac.set(N14, N14, scratch);

    jac.set(N14, N15, 0.0);

    jac.set(N14, O16, 0.0);

    scratch = screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(N14, O17, scratch);

    jac.set(N14, O18, 0.0);

    scratch = screened_rates(k_F18_to_He4_N14_derived);
    jac.set(N14, F18, scratch);

    jac.set(N14, F19, 0.0);

    jac.set(N14, Ne20, 0.0);

    jac.set(N14, Ne21, 0.0);

    jac.set(N14, Ne22, 0.0);

    jac.set(N14, Na22, 0.0);

    jac.set(N14, Na23, 0.0);

    jac.set(N14, Mg24, 0.0);

    jac.set(N14, Mg25, 0.0);

    jac.set(N14, Al26, 0.0);

    jac.set(N14, Al27, 0.0);

    jac.set(N14, Si28, 0.0);

    jac.set(N14, Si29, 0.0);

    jac.set(N14, P30, 0.0);

    jac.set(N14, P31, 0.0);

    jac.set(N14, S32, 0.0);

    jac.set(N14, S33, 0.0);

    jac.set(N14, Cl34, 0.0);

    jac.set(N14, Cl35, 0.0);

    jac.set(N14, Ar36, 0.0);

    jac.set(N14, Ar37, 0.0);

    jac.set(N14, K38, 0.0);

    jac.set(N14, K39, 0.0);

    jac.set(N14, Ca40, 0.0);

    jac.set(N14, Ca41, 0.0);

    jac.set(N14, Sc42, 0.0);

    jac.set(N14, Sc43, 0.0);

    jac.set(N14, Ti44, 0.0);

    jac.set(N14, Ti45, 0.0);

    jac.set(N14, V46, 0.0);

    jac.set(N14, V47, 0.0);

    jac.set(N14, Cr48, 0.0);

    jac.set(N14, Cr49, 0.0);

    jac.set(N14, Mn50, 0.0);

    jac.set(N14, Mn51, 0.0);

    jac.set(N14, Mn52, 0.0);

    jac.set(N14, Fe52, 0.0);

    jac.set(N14, Fe53, 0.0);

    jac.set(N14, Fe54, 0.0);

    jac.set(N14, Fe55, 0.0);

    jac.set(N14, Fe56, 0.0);

    jac.set(N14, Co55, 0.0);

    jac.set(N14, Co56, 0.0);

    jac.set(N14, Co57, 0.0);

    jac.set(N14, Ni56, 0.0);

    jac.set(N14, Ni57, 0.0);

    jac.set(N14, Ni58, 0.0);

    jac.set(N14, Cu59, 0.0);

    jac.set(N14, Zn60, 0.0);

    jac.set(N15, N, 0.0);

    scratch = screened_rates(k_p_C14_to_N15_reaclib)*Y(C14)*state.rho - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho + screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(N15, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho - screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho;
    jac.set(N15, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho;
    jac.set(N15, C12, scratch);

    jac.set(N15, C13, 0.0);

    scratch = screened_rates(k_p_C14_to_N15_reaclib)*Y(H1)*state.rho;
    jac.set(N15, C14, scratch);

    jac.set(N15, N13, 0.0);

    jac.set(N15, N14, 0.0);

    scratch = -screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho - screened_rates(k_N15_to_p_C14_derived) - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(N15, N15, scratch);

    scratch = screened_rates(k_O16_to_p_N15_derived);
    jac.set(N15, O16, scratch);

    jac.set(N15, O17, 0.0);

    scratch = screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(N15, O18, scratch);

    jac.set(N15, F18, 0.0);

    scratch = screened_rates(k_F19_to_He4_N15_derived);
    jac.set(N15, F19, scratch);

    jac.set(N15, Ne20, 0.0);

    jac.set(N15, Ne21, 0.0);

    jac.set(N15, Ne22, 0.0);

    jac.set(N15, Na22, 0.0);

    jac.set(N15, Na23, 0.0);

    jac.set(N15, Mg24, 0.0);

    jac.set(N15, Mg25, 0.0);

    jac.set(N15, Al26, 0.0);

    jac.set(N15, Al27, 0.0);

    jac.set(N15, Si28, 0.0);

    jac.set(N15, Si29, 0.0);

    jac.set(N15, P30, 0.0);

    jac.set(N15, P31, 0.0);

    jac.set(N15, S32, 0.0);

    jac.set(N15, S33, 0.0);

    jac.set(N15, Cl34, 0.0);

    jac.set(N15, Cl35, 0.0);

    jac.set(N15, Ar36, 0.0);

    jac.set(N15, Ar37, 0.0);

    jac.set(N15, K38, 0.0);

    jac.set(N15, K39, 0.0);

    jac.set(N15, Ca40, 0.0);

    jac.set(N15, Ca41, 0.0);

    jac.set(N15, Sc42, 0.0);

    jac.set(N15, Sc43, 0.0);

    jac.set(N15, Ti44, 0.0);

    jac.set(N15, Ti45, 0.0);

    jac.set(N15, V46, 0.0);

    jac.set(N15, V47, 0.0);

    jac.set(N15, Cr48, 0.0);

    jac.set(N15, Cr49, 0.0);

    jac.set(N15, Mn50, 0.0);

    jac.set(N15, Mn51, 0.0);

    jac.set(N15, Mn52, 0.0);

    jac.set(N15, Fe52, 0.0);

    jac.set(N15, Fe53, 0.0);

    jac.set(N15, Fe54, 0.0);

    jac.set(N15, Fe55, 0.0);

    jac.set(N15, Fe56, 0.0);

    jac.set(N15, Co55, 0.0);

    jac.set(N15, Co56, 0.0);

    jac.set(N15, Co57, 0.0);

    jac.set(N15, Ni56, 0.0);

    jac.set(N15, Ni57, 0.0);

    jac.set(N15, Ni58, 0.0);

    jac.set(N15, Cu59, 0.0);

    jac.set(N15, Zn60, 0.0);

    jac.set(O16, N, 0.0);

    scratch = screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*state.rho + screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho + 2.0*screened_rates(k_p_P31_to_O16_O16_derived)*Y(P31)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho + 2.0*screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(Si28)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    jac.set(O16, C13, 0.0);

    jac.set(O16, C14, 0.0);

    scratch = screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    jac.set(O16, N14, 0.0);

    scratch = screened_rates(k_p_N15_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, N15, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho - 2.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho - screened_rates(k_O16_to_He4_C12_derived) - screened_rates(k_O16_to_p_N15_derived) - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    jac.set(O16, O17, 0.0);

    jac.set(O16, O18, 0.0);

    jac.set(O16, F18, 0.0);

    scratch = screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, F19, scratch);

    scratch = screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(O16, Ne20, scratch);

    jac.set(O16, Ne21, 0.0);

    jac.set(O16, Ne22, 0.0);

    jac.set(O16, Na22, 0.0);

    jac.set(O16, Na23, 0.0);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho;
    jac.set(O16, Mg24, scratch);

    jac.set(O16, Mg25, 0.0);

    jac.set(O16, Al26, 0.0);

    scratch = screened_rates(k_p_Al27_to_C12_O16_derived)*Y(H1)*state.rho;
    jac.set(O16, Al27, scratch);

    scratch = 2.0*screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(He4)*state.rho;
    jac.set(O16, Si28, scratch);

    jac.set(O16, Si29, 0.0);

    jac.set(O16, P30, 0.0);

    scratch = 2.0*screened_rates(k_p_P31_to_O16_O16_derived)*Y(H1)*state.rho;
    jac.set(O16, P31, scratch);

    jac.set(O16, S32, 0.0);

    jac.set(O16, S33, 0.0);

    jac.set(O16, Cl34, 0.0);

    jac.set(O16, Cl35, 0.0);

    jac.set(O16, Ar36, 0.0);

    jac.set(O16, Ar37, 0.0);

    jac.set(O16, K38, 0.0);

    jac.set(O16, K39, 0.0);

    jac.set(O16, Ca40, 0.0);

    jac.set(O16, Ca41, 0.0);

    jac.set(O16, Sc42, 0.0);

    jac.set(O16, Sc43, 0.0);

    jac.set(O16, Ti44, 0.0);

    jac.set(O16, Ti45, 0.0);

    jac.set(O16, V46, 0.0);

    jac.set(O16, V47, 0.0);

    jac.set(O16, Cr48, 0.0);

    jac.set(O16, Cr49, 0.0);

    jac.set(O16, Mn50, 0.0);

    jac.set(O16, Mn51, 0.0);

    jac.set(O16, Mn52, 0.0);

    jac.set(O16, Fe52, 0.0);

    jac.set(O16, Fe53, 0.0);

    jac.set(O16, Fe54, 0.0);

    jac.set(O16, Fe55, 0.0);

    jac.set(O16, Fe56, 0.0);

    jac.set(O16, Co55, 0.0);

    jac.set(O16, Co56, 0.0);

    jac.set(O16, Co57, 0.0);

    jac.set(O16, Ni56, 0.0);

    jac.set(O16, Ni57, 0.0);

    jac.set(O16, Ni58, 0.0);

    jac.set(O16, Cu59, 0.0);

    jac.set(O16, Zn60, 0.0);

    jac.set(O17, N, 0.0);

    scratch = -screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho;
    jac.set(O17, H1, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho - screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(O17)*state.rho;
    jac.set(O17, He4, scratch);

    jac.set(O17, C12, 0.0);

    jac.set(O17, C13, 0.0);

    jac.set(O17, C14, 0.0);

    jac.set(O17, N13, 0.0);

    scratch = screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho;
    jac.set(O17, N14, scratch);

    jac.set(O17, N15, 0.0);

    jac.set(O17, O16, 0.0);

    scratch = -screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*state.rho - screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(O17, O17, scratch);

    jac.set(O17, O18, 0.0);

    scratch = screened_rates(k_F18_to_p_O17_derived);
    jac.set(O17, F18, scratch);

    jac.set(O17, F19, 0.0);

    jac.set(O17, Ne20, 0.0);

    scratch = screened_rates(k_Ne21_to_He4_O17_derived);
    jac.set(O17, Ne21, scratch);

    jac.set(O17, Ne22, 0.0);

    jac.set(O17, Na22, 0.0);

    jac.set(O17, Na23, 0.0);

    jac.set(O17, Mg24, 0.0);

    jac.set(O17, Mg25, 0.0);

    jac.set(O17, Al26, 0.0);

    jac.set(O17, Al27, 0.0);

    jac.set(O17, Si28, 0.0);

    jac.set(O17, Si29, 0.0);

    jac.set(O17, P30, 0.0);

    jac.set(O17, P31, 0.0);

    jac.set(O17, S32, 0.0);

    jac.set(O17, S33, 0.0);

    jac.set(O17, Cl34, 0.0);

    jac.set(O17, Cl35, 0.0);

    jac.set(O17, Ar36, 0.0);

    jac.set(O17, Ar37, 0.0);

    jac.set(O17, K38, 0.0);

    jac.set(O17, K39, 0.0);

    jac.set(O17, Ca40, 0.0);

    jac.set(O17, Ca41, 0.0);

    jac.set(O17, Sc42, 0.0);

    jac.set(O17, Sc43, 0.0);

    jac.set(O17, Ti44, 0.0);

    jac.set(O17, Ti45, 0.0);

    jac.set(O17, V46, 0.0);

    jac.set(O17, V47, 0.0);

    jac.set(O17, Cr48, 0.0);

    jac.set(O17, Cr49, 0.0);

    jac.set(O17, Mn50, 0.0);

    jac.set(O17, Mn51, 0.0);

    jac.set(O17, Mn52, 0.0);

    jac.set(O17, Fe52, 0.0);

    jac.set(O17, Fe53, 0.0);

    jac.set(O17, Fe54, 0.0);

    jac.set(O17, Fe55, 0.0);

    jac.set(O17, Fe56, 0.0);

    jac.set(O17, Co55, 0.0);

    jac.set(O17, Co56, 0.0);

    jac.set(O17, Co57, 0.0);

    jac.set(O17, Ni56, 0.0);

    jac.set(O17, Ni57, 0.0);

    jac.set(O17, Ni58, 0.0);

    jac.set(O17, Cu59, 0.0);

    jac.set(O17, Zn60, 0.0);

    jac.set(O18, N, 0.0);

    scratch = -screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(O18, H1, scratch);

    scratch = screened_rates(k_He4_C14_to_O18_reaclib)*Y(C14)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho - screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(O18)*state.rho;
    jac.set(O18, He4, scratch);

    jac.set(O18, C12, 0.0);

    jac.set(O18, C13, 0.0);

    scratch = screened_rates(k_He4_C14_to_O18_reaclib)*Y(He4)*state.rho;
    jac.set(O18, C14, scratch);

    jac.set(O18, N13, 0.0);

    jac.set(O18, N14, 0.0);

    scratch = screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho;
    jac.set(O18, N15, scratch);

    jac.set(O18, O16, 0.0);

    jac.set(O18, O17, 0.0);

    scratch = -screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(He4)*state.rho - screened_rates(k_O18_to_F18_weaktab) - screened_rates(k_O18_to_He4_C14_derived) - screened_rates(k_p_O18_to_F19_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(O18, O18, scratch);

    scratch = screened_rates(k_F18_to_O18_weaktab);
    jac.set(O18, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18_derived);
    jac.set(O18, F19, scratch);

    jac.set(O18, Ne20, 0.0);

    jac.set(O18, Ne21, 0.0);

    scratch = screened_rates(k_Ne22_to_He4_O18_derived);
    jac.set(O18, Ne22, scratch);

    jac.set(O18, Na22, 0.0);

    jac.set(O18, Na23, 0.0);

    jac.set(O18, Mg24, 0.0);

    jac.set(O18, Mg25, 0.0);

    jac.set(O18, Al26, 0.0);

    jac.set(O18, Al27, 0.0);

    jac.set(O18, Si28, 0.0);

    jac.set(O18, Si29, 0.0);

    jac.set(O18, P30, 0.0);

    jac.set(O18, P31, 0.0);

    jac.set(O18, S32, 0.0);

    jac.set(O18, S33, 0.0);

    jac.set(O18, Cl34, 0.0);

    jac.set(O18, Cl35, 0.0);

    jac.set(O18, Ar36, 0.0);

    jac.set(O18, Ar37, 0.0);

    jac.set(O18, K38, 0.0);

    jac.set(O18, K39, 0.0);

    jac.set(O18, Ca40, 0.0);

    jac.set(O18, Ca41, 0.0);

    jac.set(O18, Sc42, 0.0);

    jac.set(O18, Sc43, 0.0);

    jac.set(O18, Ti44, 0.0);

    jac.set(O18, Ti45, 0.0);

    jac.set(O18, V46, 0.0);

    jac.set(O18, V47, 0.0);

    jac.set(O18, Cr48, 0.0);

    jac.set(O18, Cr49, 0.0);

    jac.set(O18, Mn50, 0.0);

    jac.set(O18, Mn51, 0.0);

    jac.set(O18, Mn52, 0.0);

    jac.set(O18, Fe52, 0.0);

    jac.set(O18, Fe53, 0.0);

    jac.set(O18, Fe54, 0.0);

    jac.set(O18, Fe55, 0.0);

    jac.set(O18, Fe56, 0.0);

    jac.set(O18, Co55, 0.0);

    jac.set(O18, Co56, 0.0);

    jac.set(O18, Co57, 0.0);

    jac.set(O18, Ni56, 0.0);

    jac.set(O18, Ni57, 0.0);

    jac.set(O18, Ni58, 0.0);

    jac.set(O18, Cu59, 0.0);

    jac.set(O18, Zn60, 0.0);

    jac.set(F18, N, 0.0);

    scratch = screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho + screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = -screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho + screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho;
    jac.set(F18, He4, scratch);

    jac.set(F18, C12, 0.0);

    jac.set(F18, C13, 0.0);

    jac.set(F18, C14, 0.0);

    jac.set(F18, N13, 0.0);

    scratch = screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    jac.set(F18, N15, 0.0);

    jac.set(F18, O16, 0.0);

    scratch = screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho;
    jac.set(F18, O17, scratch);

    scratch = screened_rates(k_O18_to_F18_weaktab);
    jac.set(F18, O18, scratch);

    scratch = -screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_F18_to_O18_weaktab) - screened_rates(k_F18_to_p_O17_derived) - screened_rates(k_He4_F18_to_Na22_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho;
    jac.set(F18, F18, scratch);

    jac.set(F18, F19, 0.0);

    jac.set(F18, Ne20, 0.0);

    scratch = screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho;
    jac.set(F18, Ne21, scratch);

    jac.set(F18, Ne22, 0.0);

    scratch = screened_rates(k_Na22_to_He4_F18_derived);
    jac.set(F18, Na22, scratch);

    jac.set(F18, Na23, 0.0);

    jac.set(F18, Mg24, 0.0);

    jac.set(F18, Mg25, 0.0);

    jac.set(F18, Al26, 0.0);

    jac.set(F18, Al27, 0.0);

    jac.set(F18, Si28, 0.0);

    jac.set(F18, Si29, 0.0);

    jac.set(F18, P30, 0.0);

    jac.set(F18, P31, 0.0);

    jac.set(F18, S32, 0.0);

    jac.set(F18, S33, 0.0);

    jac.set(F18, Cl34, 0.0);

    jac.set(F18, Cl35, 0.0);

    jac.set(F18, Ar36, 0.0);

    jac.set(F18, Ar37, 0.0);

    jac.set(F18, K38, 0.0);

    jac.set(F18, K39, 0.0);

    jac.set(F18, Ca40, 0.0);

    jac.set(F18, Ca41, 0.0);

    jac.set(F18, Sc42, 0.0);

    jac.set(F18, Sc43, 0.0);

    jac.set(F18, Ti44, 0.0);

    jac.set(F18, Ti45, 0.0);

    jac.set(F18, V46, 0.0);

    jac.set(F18, V47, 0.0);

    jac.set(F18, Cr48, 0.0);

    jac.set(F18, Cr49, 0.0);

    jac.set(F18, Mn50, 0.0);

    jac.set(F18, Mn51, 0.0);

    jac.set(F18, Mn52, 0.0);

    jac.set(F18, Fe52, 0.0);

    jac.set(F18, Fe53, 0.0);

    jac.set(F18, Fe54, 0.0);

    jac.set(F18, Fe55, 0.0);

    jac.set(F18, Fe56, 0.0);

    jac.set(F18, Co55, 0.0);

    jac.set(F18, Co56, 0.0);

    jac.set(F18, Co57, 0.0);

    jac.set(F18, Ni56, 0.0);

    jac.set(F18, Ni57, 0.0);

    jac.set(F18, Ni58, 0.0);

    jac.set(F18, Cu59, 0.0);

    jac.set(F18, Zn60, 0.0);

    jac.set(F19, N, 0.0);

    scratch = -screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho + screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*state.rho + screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho;
    jac.set(F19, H1, scratch);

    scratch = -screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*state.rho - screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*state.rho + screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(F19, He4, scratch);

    jac.set(F19, C12, 0.0);

    jac.set(F19, C13, 0.0);

    jac.set(F19, C14, 0.0);

    jac.set(F19, N13, 0.0);

    jac.set(F19, N14, 0.0);

    scratch = screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*state.rho;
    jac.set(F19, N15, scratch);

    scratch = screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho;
    jac.set(F19, O16, scratch);

    jac.set(F19, O17, 0.0);

    scratch = screened_rates(k_p_O18_to_F19_reaclib)*Y(H1)*state.rho;
    jac.set(F19, O18, scratch);

    jac.set(F19, F18, 0.0);

    scratch = -screened_rates(k_F19_to_He4_N15_derived) - screened_rates(k_F19_to_p_O18_derived) - screened_rates(k_He4_F19_to_Na23_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(He4)*state.rho - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(F19, F19, scratch);

    scratch = screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(F19, Ne20, scratch);

    jac.set(F19, Ne21, 0.0);

    scratch = screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(H1)*state.rho;
    jac.set(F19, Ne22, scratch);

    jac.set(F19, Na22, 0.0);

    scratch = screened_rates(k_Na23_to_He4_F19_derived);
    jac.set(F19, Na23, scratch);

    jac.set(F19, Mg24, 0.0);

    jac.set(F19, Mg25, 0.0);

    jac.set(F19, Al26, 0.0);

    jac.set(F19, Al27, 0.0);

    jac.set(F19, Si28, 0.0);

    jac.set(F19, Si29, 0.0);

    jac.set(F19, P30, 0.0);

    jac.set(F19, P31, 0.0);

    jac.set(F19, S32, 0.0);

    jac.set(F19, S33, 0.0);

    jac.set(F19, Cl34, 0.0);

    jac.set(F19, Cl35, 0.0);

    jac.set(F19, Ar36, 0.0);

    jac.set(F19, Ar37, 0.0);

    jac.set(F19, K38, 0.0);

    jac.set(F19, K39, 0.0);

    jac.set(F19, Ca40, 0.0);

    jac.set(F19, Ca41, 0.0);

    jac.set(F19, Sc42, 0.0);

    jac.set(F19, Sc43, 0.0);

    jac.set(F19, Ti44, 0.0);

    jac.set(F19, Ti45, 0.0);

    jac.set(F19, V46, 0.0);

    jac.set(F19, V47, 0.0);

    jac.set(F19, Cr48, 0.0);

    jac.set(F19, Cr49, 0.0);

    jac.set(F19, Mn50, 0.0);

    jac.set(F19, Mn51, 0.0);

    jac.set(F19, Mn52, 0.0);

    jac.set(F19, Fe52, 0.0);

    jac.set(F19, Fe53, 0.0);

    jac.set(F19, Fe54, 0.0);

    jac.set(F19, Fe55, 0.0);

    jac.set(F19, Fe56, 0.0);

    jac.set(F19, Co55, 0.0);

    jac.set(F19, Co56, 0.0);

    jac.set(F19, Co57, 0.0);

    jac.set(F19, Ni56, 0.0);

    jac.set(F19, Ni57, 0.0);

    jac.set(F19, Ni58, 0.0);

    jac.set(F19, Cu59, 0.0);

    jac.set(F19, Zn60, 0.0);

    jac.set(Ne20, N, 0.0);

    scratch = screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho + screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho + screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(Si28)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(Ne20)*state.rho - screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(Ne20)*state.rho;
    jac.set(Ne20, C12, scratch);

    jac.set(Ne20, C13, 0.0);

    jac.set(Ne20, C14, 0.0);

    jac.set(Ne20, N13, 0.0);

    jac.set(Ne20, N14, 0.0);

    jac.set(Ne20, N15, 0.0);

    scratch = screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    jac.set(Ne20, O17, 0.0);

    jac.set(Ne20, O18, 0.0);

    jac.set(Ne20, F18, 0.0);

    scratch = screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, F19, scratch);

    scratch = -screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho - screened_rates(k_Ne20_to_He4_O16_derived) - screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(Ne20, Ne20, scratch);

    jac.set(Ne20, Ne21, 0.0);

    jac.set(Ne20, Ne22, 0.0);

    jac.set(Ne20, Na22, 0.0);

    scratch = screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, Na23, scratch);

    scratch = screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(Ne20, Mg24, scratch);

    jac.set(Ne20, Mg25, 0.0);

    jac.set(Ne20, Al26, 0.0);

    jac.set(Ne20, Al27, 0.0);

    scratch = screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*state.rho;
    jac.set(Ne20, Si28, scratch);

    jac.set(Ne20, Si29, 0.0);

    jac.set(Ne20, P30, 0.0);

    scratch = screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(H1)*state.rho;
    jac.set(Ne20, P31, scratch);

    jac.set(Ne20, S32, 0.0);

    jac.set(Ne20, S33, 0.0);

    jac.set(Ne20, Cl34, 0.0);

    jac.set(Ne20, Cl35, 0.0);

    jac.set(Ne20, Ar36, 0.0);

    jac.set(Ne20, Ar37, 0.0);

    jac.set(Ne20, K38, 0.0);

    jac.set(Ne20, K39, 0.0);

    jac.set(Ne20, Ca40, 0.0);

    jac.set(Ne20, Ca41, 0.0);

    jac.set(Ne20, Sc42, 0.0);

    jac.set(Ne20, Sc43, 0.0);

    jac.set(Ne20, Ti44, 0.0);

    jac.set(Ne20, Ti45, 0.0);

    jac.set(Ne20, V46, 0.0);

    jac.set(Ne20, V47, 0.0);

    jac.set(Ne20, Cr48, 0.0);

    jac.set(Ne20, Cr49, 0.0);

    jac.set(Ne20, Mn50, 0.0);

    jac.set(Ne20, Mn51, 0.0);

    jac.set(Ne20, Mn52, 0.0);

    jac.set(Ne20, Fe52, 0.0);

    jac.set(Ne20, Fe53, 0.0);

    jac.set(Ne20, Fe54, 0.0);

    jac.set(Ne20, Fe55, 0.0);

    jac.set(Ne20, Fe56, 0.0);

    jac.set(Ne20, Co55, 0.0);

    jac.set(Ne20, Co56, 0.0);

    jac.set(Ne20, Co57, 0.0);

    jac.set(Ne20, Ni56, 0.0);

    jac.set(Ne20, Ni57, 0.0);

    jac.set(Ne20, Ni58, 0.0);

    jac.set(Ne20, Cu59, 0.0);

    jac.set(Ne20, Zn60, 0.0);

    jac.set(Ne21, N, 0.0);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*state.rho;
    jac.set(Ne21, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(Ne21)*state.rho + screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(O17)*state.rho;
    jac.set(Ne21, He4, scratch);

    jac.set(Ne21, C12, 0.0);

    jac.set(Ne21, C13, 0.0);

    jac.set(Ne21, C14, 0.0);

    jac.set(Ne21, N13, 0.0);

    jac.set(Ne21, N14, 0.0);

    jac.set(Ne21, N15, 0.0);

    jac.set(Ne21, O16, 0.0);

    scratch = screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*state.rho;
    jac.set(Ne21, O17, scratch);

    jac.set(Ne21, O18, 0.0);

    scratch = screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho;
    jac.set(Ne21, F18, scratch);

    jac.set(Ne21, F19, 0.0);

    jac.set(Ne21, Ne20, 0.0);

    scratch = -screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(He4)*state.rho - screened_rates(k_Ne21_to_He4_O17_derived) - screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(H1)*state.rho;
    jac.set(Ne21, Ne21, scratch);

    jac.set(Ne21, Ne22, 0.0);

    scratch = screened_rates(k_Na22_to_p_Ne21_derived);
    jac.set(Ne21, Na22, scratch);

    jac.set(Ne21, Na23, 0.0);

    jac.set(Ne21, Mg24, 0.0);

    scratch = screened_rates(k_Mg25_to_He4_Ne21_derived);
    jac.set(Ne21, Mg25, scratch);

    jac.set(Ne21, Al26, 0.0);

    jac.set(Ne21, Al27, 0.0);

    jac.set(Ne21, Si28, 0.0);

    jac.set(Ne21, Si29, 0.0);

    jac.set(Ne21, P30, 0.0);

    jac.set(Ne21, P31, 0.0);

    jac.set(Ne21, S32, 0.0);

    jac.set(Ne21, S33, 0.0);

    jac.set(Ne21, Cl34, 0.0);

    jac.set(Ne21, Cl35, 0.0);

    jac.set(Ne21, Ar36, 0.0);

    jac.set(Ne21, Ar37, 0.0);

    jac.set(Ne21, K38, 0.0);

    jac.set(Ne21, K39, 0.0);

    jac.set(Ne21, Ca40, 0.0);

    jac.set(Ne21, Ca41, 0.0);

    jac.set(Ne21, Sc42, 0.0);

    jac.set(Ne21, Sc43, 0.0);

    jac.set(Ne21, Ti44, 0.0);

    jac.set(Ne21, Ti45, 0.0);

    jac.set(Ne21, V46, 0.0);

    jac.set(Ne21, V47, 0.0);

    jac.set(Ne21, Cr48, 0.0);

    jac.set(Ne21, Cr49, 0.0);

    jac.set(Ne21, Mn50, 0.0);

    jac.set(Ne21, Mn51, 0.0);

    jac.set(Ne21, Mn52, 0.0);

    jac.set(Ne21, Fe52, 0.0);

    jac.set(Ne21, Fe53, 0.0);

    jac.set(Ne21, Fe54, 0.0);

    jac.set(Ne21, Fe55, 0.0);

    jac.set(Ne21, Fe56, 0.0);

    jac.set(Ne21, Co55, 0.0);

    jac.set(Ne21, Co56, 0.0);

    jac.set(Ne21, Co57, 0.0);

    jac.set(Ne21, Ni56, 0.0);

    jac.set(Ne21, Ni57, 0.0);

    jac.set(Ne21, Ni58, 0.0);

    jac.set(Ne21, Cu59, 0.0);

    jac.set(Ne21, Zn60, 0.0);

    jac.set(Ne22, N, 0.0);

    scratch = -screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(Ne22)*state.rho - screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(Ne22)*state.rho;
    jac.set(Ne22, H1, scratch);

    scratch = screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(F19)*state.rho + screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(O18)*state.rho;
    jac.set(Ne22, He4, scratch);

    jac.set(Ne22, C12, 0.0);

    jac.set(Ne22, C13, 0.0);

    jac.set(Ne22, C14, 0.0);

    jac.set(Ne22, N13, 0.0);

    jac.set(Ne22, N14, 0.0);

    jac.set(Ne22, N15, 0.0);

    jac.set(Ne22, O16, 0.0);

    jac.set(Ne22, O17, 0.0);

    scratch = screened_rates(k_He4_O18_to_Ne22_reaclib)*Y(He4)*state.rho;
    jac.set(Ne22, O18, scratch);

    jac.set(Ne22, F18, 0.0);

    scratch = screened_rates(k_He4_F19_to_p_Ne22_reaclib)*Y(He4)*state.rho;
    jac.set(Ne22, F19, scratch);

    jac.set(Ne22, Ne20, 0.0);

    jac.set(Ne22, Ne21, 0.0);

    scratch = -screened_rates(k_Ne22_to_He4_O18_derived) - screened_rates(k_Ne22_to_Na22_weaktab) - screened_rates(k_p_Ne22_to_He4_F19_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(H1)*state.rho;
    jac.set(Ne22, Ne22, scratch);

    scratch = screened_rates(k_Na22_to_Ne22_weaktab);
    jac.set(Ne22, Na22, scratch);

    scratch = screened_rates(k_Na23_to_p_Ne22_derived);
    jac.set(Ne22, Na23, scratch);

    jac.set(Ne22, Mg24, 0.0);

    jac.set(Ne22, Mg25, 0.0);

    jac.set(Ne22, Al26, 0.0);

    jac.set(Ne22, Al27, 0.0);

    jac.set(Ne22, Si28, 0.0);

    jac.set(Ne22, Si29, 0.0);

    jac.set(Ne22, P30, 0.0);

    jac.set(Ne22, P31, 0.0);

    jac.set(Ne22, S32, 0.0);

    jac.set(Ne22, S33, 0.0);

    jac.set(Ne22, Cl34, 0.0);

    jac.set(Ne22, Cl35, 0.0);

    jac.set(Ne22, Ar36, 0.0);

    jac.set(Ne22, Ar37, 0.0);

    jac.set(Ne22, K38, 0.0);

    jac.set(Ne22, K39, 0.0);

    jac.set(Ne22, Ca40, 0.0);

    jac.set(Ne22, Ca41, 0.0);

    jac.set(Ne22, Sc42, 0.0);

    jac.set(Ne22, Sc43, 0.0);

    jac.set(Ne22, Ti44, 0.0);

    jac.set(Ne22, Ti45, 0.0);

    jac.set(Ne22, V46, 0.0);

    jac.set(Ne22, V47, 0.0);

    jac.set(Ne22, Cr48, 0.0);

    jac.set(Ne22, Cr49, 0.0);

    jac.set(Ne22, Mn50, 0.0);

    jac.set(Ne22, Mn51, 0.0);

    jac.set(Ne22, Mn52, 0.0);

    jac.set(Ne22, Fe52, 0.0);

    jac.set(Ne22, Fe53, 0.0);

    jac.set(Ne22, Fe54, 0.0);

    jac.set(Ne22, Fe55, 0.0);

    jac.set(Ne22, Fe56, 0.0);

    jac.set(Ne22, Co55, 0.0);

    jac.set(Ne22, Co56, 0.0);

    jac.set(Ne22, Co57, 0.0);

    jac.set(Ne22, Ni56, 0.0);

    jac.set(Ne22, Ni57, 0.0);

    jac.set(Ne22, Ni58, 0.0);

    jac.set(Ne22, Cu59, 0.0);

    jac.set(Ne22, Zn60, 0.0);

    jac.set(Na22, N, 0.0);

    scratch = screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*state.rho + screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*state.rho;
    jac.set(Na22, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(Na22)*state.rho - screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(Na22)*state.rho;
    jac.set(Na22, He4, scratch);

    jac.set(Na22, C12, 0.0);

    jac.set(Na22, C13, 0.0);

    jac.set(Na22, C14, 0.0);

    jac.set(Na22, N13, 0.0);

    jac.set(Na22, N14, 0.0);

    jac.set(Na22, N15, 0.0);

    jac.set(Na22, O16, 0.0);

    jac.set(Na22, O17, 0.0);

    jac.set(Na22, O18, 0.0);

    scratch = screened_rates(k_He4_F18_to_Na22_reaclib)*Y(He4)*state.rho;
    jac.set(Na22, F18, scratch);

    jac.set(Na22, F19, 0.0);

    jac.set(Na22, Ne20, 0.0);

    scratch = screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(H1)*state.rho;
    jac.set(Na22, Ne21, scratch);

    scratch = screened_rates(k_Ne22_to_Na22_weaktab);
    jac.set(Na22, Ne22, scratch);

    scratch = -screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*state.rho - screened_rates(k_Na22_to_He4_F18_derived) - screened_rates(k_Na22_to_Ne22_weaktab) - screened_rates(k_Na22_to_p_Ne21_derived);
    jac.set(Na22, Na22, scratch);

    jac.set(Na22, Na23, 0.0);

    jac.set(Na22, Mg24, 0.0);

    scratch = screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(H1)*state.rho;
    jac.set(Na22, Mg25, scratch);

    scratch = screened_rates(k_Al26_to_He4_Na22_derived);
    jac.set(Na22, Al26, scratch);

    jac.set(Na22, Al27, 0.0);

    jac.set(Na22, Si28, 0.0);

    jac.set(Na22, Si29, 0.0);

    jac.set(Na22, P30, 0.0);

    jac.set(Na22, P31, 0.0);

    jac.set(Na22, S32, 0.0);

    jac.set(Na22, S33, 0.0);

    jac.set(Na22, Cl34, 0.0);

    jac.set(Na22, Cl35, 0.0);

    jac.set(Na22, Ar36, 0.0);

    jac.set(Na22, Ar37, 0.0);

    jac.set(Na22, K38, 0.0);

    jac.set(Na22, K39, 0.0);

    jac.set(Na22, Ca40, 0.0);

    jac.set(Na22, Ca41, 0.0);

    jac.set(Na22, Sc42, 0.0);

    jac.set(Na22, Sc43, 0.0);

    jac.set(Na22, Ti44, 0.0);

    jac.set(Na22, Ti45, 0.0);

    jac.set(Na22, V46, 0.0);

    jac.set(Na22, V47, 0.0);

    jac.set(Na22, Cr48, 0.0);

    jac.set(Na22, Cr49, 0.0);

    jac.set(Na22, Mn50, 0.0);

    jac.set(Na22, Mn51, 0.0);

    jac.set(Na22, Mn52, 0.0);

    jac.set(Na22, Fe52, 0.0);

    jac.set(Na22, Fe53, 0.0);

    jac.set(Na22, Fe54, 0.0);

    jac.set(Na22, Fe55, 0.0);

    jac.set(Na22, Fe56, 0.0);

    jac.set(Na22, Co55, 0.0);

    jac.set(Na22, Co56, 0.0);

    jac.set(Na22, Co57, 0.0);

    jac.set(Na22, Ni56, 0.0);

    jac.set(Na22, Ni57, 0.0);

    jac.set(Na22, Ni58, 0.0);

    jac.set(Na22, Cu59, 0.0);

    jac.set(Na22, Zn60, 0.0);

    jac.set(Na23, N, 0.0);

    scratch = -screened_rates(k_p_Na23_to_C12_C12_derived)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*state.rho + screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(Ne22)*state.rho;
    jac.set(Na23, H1, scratch);

    scratch = screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*state.rho - screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(Na23)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho;
    jac.set(Na23, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    jac.set(Na23, C13, 0.0);

    jac.set(Na23, C14, 0.0);

    jac.set(Na23, N13, 0.0);

    jac.set(Na23, N14, 0.0);

    jac.set(Na23, N15, 0.0);

    jac.set(Na23, O16, 0.0);

    jac.set(Na23, O17, 0.0);

    jac.set(Na23, O18, 0.0);

    jac.set(Na23, F18, 0.0);

    scratch = screened_rates(k_He4_F19_to_Na23_reaclib)*Y(He4)*state.rho;
    jac.set(Na23, F19, scratch);

    scratch = screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho;
    jac.set(Na23, Ne20, scratch);

    jac.set(Na23, Ne21, 0.0);

    scratch = screened_rates(k_p_Ne22_to_Na23_reaclib)*Y(H1)*state.rho;
    jac.set(Na23, Ne22, scratch);

    jac.set(Na23, Na22, 0.0);

    scratch = -screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(He4)*state.rho - screened_rates(k_Na23_to_He4_F19_derived) - screened_rates(k_Na23_to_p_Ne22_derived) - screened_rates(k_p_Na23_to_C12_C12_derived)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Na23, Na23, scratch);

    scratch = screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(Na23, Mg24, scratch);

    jac.set(Na23, Mg25, 0.0);

    jac.set(Na23, Al26, 0.0);

    scratch = screened_rates(k_Al27_to_He4_Na23_derived);
    jac.set(Na23, Al27, scratch);

    jac.set(Na23, Si28, 0.0);

    jac.set(Na23, Si29, 0.0);

    jac.set(Na23, P30, 0.0);

    jac.set(Na23, P31, 0.0);

    jac.set(Na23, S32, 0.0);

    jac.set(Na23, S33, 0.0);

    jac.set(Na23, Cl34, 0.0);

    jac.set(Na23, Cl35, 0.0);

    jac.set(Na23, Ar36, 0.0);

    jac.set(Na23, Ar37, 0.0);

    jac.set(Na23, K38, 0.0);

    jac.set(Na23, K39, 0.0);

    jac.set(Na23, Ca40, 0.0);

    jac.set(Na23, Ca41, 0.0);

    jac.set(Na23, Sc42, 0.0);

    jac.set(Na23, Sc43, 0.0);

    jac.set(Na23, Ti44, 0.0);

    jac.set(Na23, Ti45, 0.0);

    jac.set(Na23, V46, 0.0);

    jac.set(Na23, V47, 0.0);

    jac.set(Na23, Cr48, 0.0);

    jac.set(Na23, Cr49, 0.0);

    jac.set(Na23, Mn50, 0.0);

    jac.set(Na23, Mn51, 0.0);

    jac.set(Na23, Mn52, 0.0);

    jac.set(Na23, Fe52, 0.0);

    jac.set(Na23, Fe53, 0.0);

    jac.set(Na23, Fe54, 0.0);

    jac.set(Na23, Fe55, 0.0);

    jac.set(Na23, Fe56, 0.0);

    jac.set(Na23, Co55, 0.0);

    jac.set(Na23, Co56, 0.0);

    jac.set(Na23, Co57, 0.0);

    jac.set(Na23, Ni56, 0.0);

    jac.set(Na23, Ni57, 0.0);

    jac.set(Na23, Ni58, 0.0);

    jac.set(Na23, Cu59, 0.0);

    jac.set(Na23, Zn60, 0.0);

    jac.set(Mg24, N, 0.0);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho + screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    jac.set(Mg24, C13, 0.0);

    jac.set(Mg24, C14, 0.0);

    jac.set(Mg24, N13, 0.0);

    jac.set(Mg24, N14, 0.0);

    jac.set(Mg24, N15, 0.0);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    jac.set(Mg24, O17, 0.0);

    jac.set(Mg24, O18, 0.0);

    jac.set(Mg24, F18, 0.0);

    jac.set(Mg24, F19, 0.0);

    scratch = screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    jac.set(Mg24, Ne21, 0.0);

    jac.set(Mg24, Ne22, 0.0);

    jac.set(Mg24, Na22, 0.0);

    scratch = screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Mg24, Na23, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho - screened_rates(k_Mg24_to_He4_Ne20_derived) - screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(Mg24, Mg24, scratch);

    jac.set(Mg24, Mg25, 0.0);

    jac.set(Mg24, Al26, 0.0);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    scratch = screened_rates(k_Si28_to_He4_Mg24_derived);
    jac.set(Mg24, Si28, scratch);

    jac.set(Mg24, Si29, 0.0);

    jac.set(Mg24, P30, 0.0);

    jac.set(Mg24, P31, 0.0);

    jac.set(Mg24, S32, 0.0);

    jac.set(Mg24, S33, 0.0);

    jac.set(Mg24, Cl34, 0.0);

    jac.set(Mg24, Cl35, 0.0);

    jac.set(Mg24, Ar36, 0.0);

    jac.set(Mg24, Ar37, 0.0);

    jac.set(Mg24, K38, 0.0);

    jac.set(Mg24, K39, 0.0);

    jac.set(Mg24, Ca40, 0.0);

    jac.set(Mg24, Ca41, 0.0);

    jac.set(Mg24, Sc42, 0.0);

    jac.set(Mg24, Sc43, 0.0);

    jac.set(Mg24, Ti44, 0.0);

    jac.set(Mg24, Ti45, 0.0);

    jac.set(Mg24, V46, 0.0);

    jac.set(Mg24, V47, 0.0);

    jac.set(Mg24, Cr48, 0.0);

    jac.set(Mg24, Cr49, 0.0);

    jac.set(Mg24, Mn50, 0.0);

    jac.set(Mg24, Mn51, 0.0);

    jac.set(Mg24, Mn52, 0.0);

    jac.set(Mg24, Fe52, 0.0);

    jac.set(Mg24, Fe53, 0.0);

    jac.set(Mg24, Fe54, 0.0);

    jac.set(Mg24, Fe55, 0.0);

    jac.set(Mg24, Fe56, 0.0);

    jac.set(Mg24, Co55, 0.0);

    jac.set(Mg24, Co56, 0.0);

    jac.set(Mg24, Co57, 0.0);

    jac.set(Mg24, Ni56, 0.0);

    jac.set(Mg24, Ni57, 0.0);

    jac.set(Mg24, Ni58, 0.0);

    jac.set(Mg24, Cu59, 0.0);

    jac.set(Mg24, Zn60, 0.0);

    jac.set(Mg25, N, 0.0);

    scratch = -screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(Mg25)*state.rho - screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(Mg25)*state.rho;
    jac.set(Mg25, H1, scratch);

    scratch = -screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(Mg25)*state.rho + screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(Na22)*state.rho + screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(Ne21)*state.rho;
    jac.set(Mg25, He4, scratch);

    jac.set(Mg25, C12, 0.0);

    jac.set(Mg25, C13, 0.0);

    jac.set(Mg25, C14, 0.0);

    jac.set(Mg25, N13, 0.0);

    jac.set(Mg25, N14, 0.0);

    jac.set(Mg25, N15, 0.0);

    jac.set(Mg25, O16, 0.0);

    jac.set(Mg25, O17, 0.0);

    jac.set(Mg25, O18, 0.0);

    jac.set(Mg25, F18, 0.0);

    jac.set(Mg25, F19, 0.0);

    jac.set(Mg25, Ne20, 0.0);

    scratch = screened_rates(k_He4_Ne21_to_Mg25_reaclib)*Y(He4)*state.rho;
    jac.set(Mg25, Ne21, scratch);

    jac.set(Mg25, Ne22, 0.0);

    scratch = screened_rates(k_He4_Na22_to_p_Mg25_reaclib)*Y(He4)*state.rho;
    jac.set(Mg25, Na22, scratch);

    jac.set(Mg25, Na23, 0.0);

    jac.set(Mg25, Mg24, 0.0);

    scratch = -screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(He4)*state.rho - screened_rates(k_Mg25_to_He4_Ne21_derived) - screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Mg25_to_He4_Na22_derived)*Y(H1)*state.rho;
    jac.set(Mg25, Mg25, scratch);

    scratch = screened_rates(k_Al26_to_p_Mg25_derived);
    jac.set(Mg25, Al26, scratch);

    jac.set(Mg25, Al27, 0.0);

    jac.set(Mg25, Si28, 0.0);

    scratch = screened_rates(k_Si29_to_He4_Mg25_derived);
    jac.set(Mg25, Si29, scratch);

    jac.set(Mg25, P30, 0.0);

    jac.set(Mg25, P31, 0.0);

    jac.set(Mg25, S32, 0.0);

    jac.set(Mg25, S33, 0.0);

    jac.set(Mg25, Cl34, 0.0);

    jac.set(Mg25, Cl35, 0.0);

    jac.set(Mg25, Ar36, 0.0);

    jac.set(Mg25, Ar37, 0.0);

    jac.set(Mg25, K38, 0.0);

    jac.set(Mg25, K39, 0.0);

    jac.set(Mg25, Ca40, 0.0);

    jac.set(Mg25, Ca41, 0.0);

    jac.set(Mg25, Sc42, 0.0);

    jac.set(Mg25, Sc43, 0.0);

    jac.set(Mg25, Ti44, 0.0);

    jac.set(Mg25, Ti45, 0.0);

    jac.set(Mg25, V46, 0.0);

    jac.set(Mg25, V47, 0.0);

    jac.set(Mg25, Cr48, 0.0);

    jac.set(Mg25, Cr49, 0.0);

    jac.set(Mg25, Mn50, 0.0);

    jac.set(Mg25, Mn51, 0.0);

    jac.set(Mg25, Mn52, 0.0);

    jac.set(Mg25, Fe52, 0.0);

    jac.set(Mg25, Fe53, 0.0);

    jac.set(Mg25, Fe54, 0.0);

    jac.set(Mg25, Fe55, 0.0);

    jac.set(Mg25, Fe56, 0.0);

    jac.set(Mg25, Co55, 0.0);

    jac.set(Mg25, Co56, 0.0);

    jac.set(Mg25, Co57, 0.0);

    jac.set(Mg25, Ni56, 0.0);

    jac.set(Mg25, Ni57, 0.0);

    jac.set(Mg25, Ni58, 0.0);

    jac.set(Mg25, Cu59, 0.0);

    jac.set(Mg25, Zn60, 0.0);

    jac.set(Al26, N, 0.0);

    scratch = screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(Mg25)*state.rho + screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*state.rho;
    jac.set(Al26, H1, scratch);

    scratch = -screened_rates(k_He4_Al26_to_P30_reaclib)*Y(Al26)*state.rho - screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*state.rho + screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(Na22)*state.rho;
    jac.set(Al26, He4, scratch);

    jac.set(Al26, C12, 0.0);

    jac.set(Al26, C13, 0.0);

    jac.set(Al26, C14, 0.0);

    jac.set(Al26, N13, 0.0);

    jac.set(Al26, N14, 0.0);

    jac.set(Al26, N15, 0.0);

    jac.set(Al26, O16, 0.0);

    jac.set(Al26, O17, 0.0);

    jac.set(Al26, O18, 0.0);

    jac.set(Al26, F18, 0.0);

    jac.set(Al26, F19, 0.0);

    jac.set(Al26, Ne20, 0.0);

    jac.set(Al26, Ne21, 0.0);

    jac.set(Al26, Ne22, 0.0);

    scratch = screened_rates(k_He4_Na22_to_Al26_reaclib)*Y(He4)*state.rho;
    jac.set(Al26, Na22, scratch);

    jac.set(Al26, Na23, 0.0);

    jac.set(Al26, Mg24, 0.0);

    scratch = screened_rates(k_p_Mg25_to_Al26_reaclib)*Y(H1)*state.rho;
    jac.set(Al26, Mg25, scratch);

    scratch = -screened_rates(k_Al26_to_He4_Na22_derived) - screened_rates(k_Al26_to_p_Mg25_derived) - screened_rates(k_He4_Al26_to_P30_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(He4)*state.rho;
    jac.set(Al26, Al26, scratch);

    jac.set(Al26, Al27, 0.0);

    jac.set(Al26, Si28, 0.0);

    scratch = screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(H1)*state.rho;
    jac.set(Al26, Si29, scratch);

    scratch = screened_rates(k_P30_to_He4_Al26_derived);
    jac.set(Al26, P30, scratch);

    jac.set(Al26, P31, 0.0);

    jac.set(Al26, S32, 0.0);

    jac.set(Al26, S33, 0.0);

    jac.set(Al26, Cl34, 0.0);

    jac.set(Al26, Cl35, 0.0);

    jac.set(Al26, Ar36, 0.0);

    jac.set(Al26, Ar37, 0.0);

    jac.set(Al26, K38, 0.0);

    jac.set(Al26, K39, 0.0);

    jac.set(Al26, Ca40, 0.0);

    jac.set(Al26, Ca41, 0.0);

    jac.set(Al26, Sc42, 0.0);

    jac.set(Al26, Sc43, 0.0);

    jac.set(Al26, Ti44, 0.0);

    jac.set(Al26, Ti45, 0.0);

    jac.set(Al26, V46, 0.0);

    jac.set(Al26, V47, 0.0);

    jac.set(Al26, Cr48, 0.0);

    jac.set(Al26, Cr49, 0.0);

    jac.set(Al26, Mn50, 0.0);

    jac.set(Al26, Mn51, 0.0);

    jac.set(Al26, Mn52, 0.0);

    jac.set(Al26, Fe52, 0.0);

    jac.set(Al26, Fe53, 0.0);

    jac.set(Al26, Fe54, 0.0);

    jac.set(Al26, Fe55, 0.0);

    jac.set(Al26, Fe56, 0.0);

    jac.set(Al26, Co55, 0.0);

    jac.set(Al26, Co56, 0.0);

    jac.set(Al26, Co57, 0.0);

    jac.set(Al26, Ni56, 0.0);

    jac.set(Al26, Ni57, 0.0);

    jac.set(Al26, Ni58, 0.0);

    jac.set(Al26, Cu59, 0.0);

    jac.set(Al26, Zn60, 0.0);

    jac.set(Al27, N, 0.0);

    scratch = -screened_rates(k_p_Al27_to_C12_O16_derived)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = -screened_rates(k_He4_Al27_to_P31_reaclib)*Y(Al27)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(Na23)*state.rho;
    jac.set(Al27, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho;
    jac.set(Al27, C12, scratch);

    jac.set(Al27, C13, 0.0);

    jac.set(Al27, C14, 0.0);

    jac.set(Al27, N13, 0.0);

    jac.set(Al27, N14, 0.0);

    jac.set(Al27, N15, 0.0);

    scratch = screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho;
    jac.set(Al27, O16, scratch);

    jac.set(Al27, O17, 0.0);

    jac.set(Al27, O18, 0.0);

    jac.set(Al27, F18, 0.0);

    jac.set(Al27, F19, 0.0);

    jac.set(Al27, Ne20, 0.0);

    jac.set(Al27, Ne21, 0.0);

    jac.set(Al27, Ne22, 0.0);

    jac.set(Al27, Na22, 0.0);

    scratch = screened_rates(k_He4_Na23_to_Al27_reaclib)*Y(He4)*state.rho;
    jac.set(Al27, Na23, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    jac.set(Al27, Mg25, 0.0);

    jac.set(Al27, Al26, 0.0);

    scratch = -screened_rates(k_Al27_to_He4_Na23_derived) - screened_rates(k_He4_Al27_to_P31_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Al27_to_C12_O16_derived)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Al27, Al27, scratch);

    scratch = screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(Al27, Si28, scratch);

    jac.set(Al27, Si29, 0.0);

    jac.set(Al27, P30, 0.0);

    scratch = screened_rates(k_P31_to_He4_Al27_derived);
    jac.set(Al27, P31, scratch);

    jac.set(Al27, S32, 0.0);

    jac.set(Al27, S33, 0.0);

    jac.set(Al27, Cl34, 0.0);

    jac.set(Al27, Cl35, 0.0);

    jac.set(Al27, Ar36, 0.0);

    jac.set(Al27, Ar37, 0.0);

    jac.set(Al27, K38, 0.0);

    jac.set(Al27, K39, 0.0);

    jac.set(Al27, Ca40, 0.0);

    jac.set(Al27, Ca41, 0.0);

    jac.set(Al27, Sc42, 0.0);

    jac.set(Al27, Sc43, 0.0);

    jac.set(Al27, Ti44, 0.0);

    jac.set(Al27, Ti45, 0.0);

    jac.set(Al27, V46, 0.0);

    jac.set(Al27, V47, 0.0);

    jac.set(Al27, Cr48, 0.0);

    jac.set(Al27, Cr49, 0.0);

    jac.set(Al27, Mn50, 0.0);

    jac.set(Al27, Mn51, 0.0);

    jac.set(Al27, Mn52, 0.0);

    jac.set(Al27, Fe52, 0.0);

    jac.set(Al27, Fe53, 0.0);

    jac.set(Al27, Fe54, 0.0);

    jac.set(Al27, Fe55, 0.0);

    jac.set(Al27, Fe56, 0.0);

    jac.set(Al27, Co55, 0.0);

    jac.set(Al27, Co56, 0.0);

    jac.set(Al27, Co57, 0.0);

    jac.set(Al27, Ni56, 0.0);

    jac.set(Al27, Ni57, 0.0);

    jac.set(Al27, Ni58, 0.0);

    jac.set(Al27, Cu59, 0.0);

    jac.set(Al27, Zn60, 0.0);

    jac.set(Si28, N, 0.0);

    scratch = screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*state.rho + screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_S32_reaclib)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(Ne20)*state.rho;
    jac.set(Si28, C12, scratch);

    jac.set(Si28, C13, 0.0);

    jac.set(Si28, C14, 0.0);

    jac.set(Si28, N13, 0.0);

    jac.set(Si28, N14, 0.0);

    jac.set(Si28, N15, 0.0);

    scratch = 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    jac.set(Si28, O17, 0.0);

    jac.set(Si28, O18, 0.0);

    jac.set(Si28, F18, 0.0);

    jac.set(Si28, F19, 0.0);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28_reaclib)*Y(C12)*state.rho;
    jac.set(Si28, Ne20, scratch);

    jac.set(Si28, Ne21, 0.0);

    jac.set(Si28, Ne22, 0.0);

    jac.set(Si28, Na22, 0.0);

    jac.set(Si28, Na23, 0.0);

    scratch = screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    jac.set(Si28, Mg25, 0.0);

    jac.set(Si28, Al26, 0.0);

    scratch = screened_rates(k_p_Al27_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_C12_Ne20_derived)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_O16_O16_derived)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho - screened_rates(k_Si28_to_He4_Mg24_derived) - screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(Si28, Si28, scratch);

    jac.set(Si28, Si29, 0.0);

    jac.set(Si28, P30, 0.0);

    scratch = screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    scratch = screened_rates(k_S32_to_He4_Si28_derived);
    jac.set(Si28, S32, scratch);

    jac.set(Si28, S33, 0.0);

    jac.set(Si28, Cl34, 0.0);

    jac.set(Si28, Cl35, 0.0);

    jac.set(Si28, Ar36, 0.0);

    jac.set(Si28, Ar37, 0.0);

    jac.set(Si28, K38, 0.0);

    jac.set(Si28, K39, 0.0);

    jac.set(Si28, Ca40, 0.0);

    jac.set(Si28, Ca41, 0.0);

    jac.set(Si28, Sc42, 0.0);

    jac.set(Si28, Sc43, 0.0);

    jac.set(Si28, Ti44, 0.0);

    jac.set(Si28, Ti45, 0.0);

    jac.set(Si28, V46, 0.0);

    jac.set(Si28, V47, 0.0);

    jac.set(Si28, Cr48, 0.0);

    jac.set(Si28, Cr49, 0.0);

    jac.set(Si28, Mn50, 0.0);

    jac.set(Si28, Mn51, 0.0);

    jac.set(Si28, Mn52, 0.0);

    jac.set(Si28, Fe52, 0.0);

    jac.set(Si28, Fe53, 0.0);

    jac.set(Si28, Fe54, 0.0);

    jac.set(Si28, Fe55, 0.0);

    jac.set(Si28, Fe56, 0.0);

    jac.set(Si28, Co55, 0.0);

    jac.set(Si28, Co56, 0.0);

    jac.set(Si28, Co57, 0.0);

    jac.set(Si28, Ni56, 0.0);

    jac.set(Si28, Ni57, 0.0);

    jac.set(Si28, Ni58, 0.0);

    jac.set(Si28, Cu59, 0.0);

    jac.set(Si28, Zn60, 0.0);

    jac.set(Si29, N, 0.0);

    scratch = -screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(Si29)*state.rho - screened_rates(k_p_Si29_to_P30_reaclib)*Y(Si29)*state.rho;
    jac.set(Si29, H1, scratch);

    scratch = screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(Al26)*state.rho + screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(Mg25)*state.rho - screened_rates(k_He4_Si29_to_S33_reaclib)*Y(Si29)*state.rho;
    jac.set(Si29, He4, scratch);

    jac.set(Si29, C12, 0.0);

    jac.set(Si29, C13, 0.0);

    jac.set(Si29, C14, 0.0);

    jac.set(Si29, N13, 0.0);

    jac.set(Si29, N14, 0.0);

    jac.set(Si29, N15, 0.0);

    jac.set(Si29, O16, 0.0);

    jac.set(Si29, O17, 0.0);

    jac.set(Si29, O18, 0.0);

    jac.set(Si29, F18, 0.0);

    jac.set(Si29, F19, 0.0);

    jac.set(Si29, Ne20, 0.0);

    jac.set(Si29, Ne21, 0.0);

    jac.set(Si29, Ne22, 0.0);

    jac.set(Si29, Na22, 0.0);

    jac.set(Si29, Na23, 0.0);

    jac.set(Si29, Mg24, 0.0);

    scratch = screened_rates(k_He4_Mg25_to_Si29_reaclib)*Y(He4)*state.rho;
    jac.set(Si29, Mg25, scratch);

    scratch = screened_rates(k_He4_Al26_to_p_Si29_reaclib)*Y(He4)*state.rho;
    jac.set(Si29, Al26, scratch);

    jac.set(Si29, Al27, 0.0);

    jac.set(Si29, Si28, 0.0);

    scratch = -screened_rates(k_He4_Si29_to_S33_reaclib)*Y(He4)*state.rho - screened_rates(k_Si29_to_He4_Mg25_derived) - screened_rates(k_p_Si29_to_He4_Al26_derived)*Y(H1)*state.rho - screened_rates(k_p_Si29_to_P30_reaclib)*Y(H1)*state.rho;
    jac.set(Si29, Si29, scratch);

    scratch = screened_rates(k_P30_to_p_Si29_derived);
    jac.set(Si29, P30, scratch);

    jac.set(Si29, P31, 0.0);

    jac.set(Si29, S32, 0.0);

    scratch = screened_rates(k_S33_to_He4_Si29_derived);
    jac.set(Si29, S33, scratch);

    jac.set(Si29, Cl34, 0.0);

    jac.set(Si29, Cl35, 0.0);

    jac.set(Si29, Ar36, 0.0);

    jac.set(Si29, Ar37, 0.0);

    jac.set(Si29, K38, 0.0);

    jac.set(Si29, K39, 0.0);

    jac.set(Si29, Ca40, 0.0);

    jac.set(Si29, Ca41, 0.0);

    jac.set(Si29, Sc42, 0.0);

    jac.set(Si29, Sc43, 0.0);

    jac.set(Si29, Ti44, 0.0);

    jac.set(Si29, Ti45, 0.0);

    jac.set(Si29, V46, 0.0);

    jac.set(Si29, V47, 0.0);

    jac.set(Si29, Cr48, 0.0);

    jac.set(Si29, Cr49, 0.0);

    jac.set(Si29, Mn50, 0.0);

    jac.set(Si29, Mn51, 0.0);

    jac.set(Si29, Mn52, 0.0);

    jac.set(Si29, Fe52, 0.0);

    jac.set(Si29, Fe53, 0.0);

    jac.set(Si29, Fe54, 0.0);

    jac.set(Si29, Fe55, 0.0);

    jac.set(Si29, Fe56, 0.0);

    jac.set(Si29, Co55, 0.0);

    jac.set(Si29, Co56, 0.0);

    jac.set(Si29, Co57, 0.0);

    jac.set(Si29, Ni56, 0.0);

    jac.set(Si29, Ni57, 0.0);

    jac.set(Si29, Ni58, 0.0);

    jac.set(Si29, Cu59, 0.0);

    jac.set(Si29, Zn60, 0.0);

    jac.set(P30, N, 0.0);

    scratch = screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*state.rho + screened_rates(k_p_Si29_to_P30_reaclib)*Y(Si29)*state.rho;
    jac.set(P30, H1, scratch);

    scratch = screened_rates(k_He4_Al26_to_P30_reaclib)*Y(Al26)*state.rho - screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(P30)*state.rho - screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(P30)*state.rho;
    jac.set(P30, He4, scratch);

    jac.set(P30, C12, 0.0);

    jac.set(P30, C13, 0.0);

    jac.set(P30, C14, 0.0);

    jac.set(P30, N13, 0.0);

    jac.set(P30, N14, 0.0);

    jac.set(P30, N15, 0.0);

    jac.set(P30, O16, 0.0);

    jac.set(P30, O17, 0.0);

    jac.set(P30, O18, 0.0);

    jac.set(P30, F18, 0.0);

    jac.set(P30, F19, 0.0);

    jac.set(P30, Ne20, 0.0);

    jac.set(P30, Ne21, 0.0);

    jac.set(P30, Ne22, 0.0);

    jac.set(P30, Na22, 0.0);

    jac.set(P30, Na23, 0.0);

    jac.set(P30, Mg24, 0.0);

    jac.set(P30, Mg25, 0.0);

    scratch = screened_rates(k_He4_Al26_to_P30_reaclib)*Y(He4)*state.rho;
    jac.set(P30, Al26, scratch);

    jac.set(P30, Al27, 0.0);

    jac.set(P30, Si28, 0.0);

    scratch = screened_rates(k_p_Si29_to_P30_reaclib)*Y(H1)*state.rho;
    jac.set(P30, Si29, scratch);

    scratch = -screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*state.rho - screened_rates(k_P30_to_He4_Al26_derived) - screened_rates(k_P30_to_p_Si29_derived);
    jac.set(P30, P30, scratch);

    jac.set(P30, P31, 0.0);

    jac.set(P30, S32, 0.0);

    scratch = screened_rates(k_p_S33_to_He4_P30_derived)*Y(H1)*state.rho;
    jac.set(P30, S33, scratch);

    scratch = screened_rates(k_Cl34_to_He4_P30_derived);
    jac.set(P30, Cl34, scratch);

    jac.set(P30, Cl35, 0.0);

    jac.set(P30, Ar36, 0.0);

    jac.set(P30, Ar37, 0.0);

    jac.set(P30, K38, 0.0);

    jac.set(P30, K39, 0.0);

    jac.set(P30, Ca40, 0.0);

    jac.set(P30, Ca41, 0.0);

    jac.set(P30, Sc42, 0.0);

    jac.set(P30, Sc43, 0.0);

    jac.set(P30, Ti44, 0.0);

    jac.set(P30, Ti45, 0.0);

    jac.set(P30, V46, 0.0);

    jac.set(P30, V47, 0.0);

    jac.set(P30, Cr48, 0.0);

    jac.set(P30, Cr49, 0.0);

    jac.set(P30, Mn50, 0.0);

    jac.set(P30, Mn51, 0.0);

    jac.set(P30, Mn52, 0.0);

    jac.set(P30, Fe52, 0.0);

    jac.set(P30, Fe53, 0.0);

    jac.set(P30, Fe54, 0.0);

    jac.set(P30, Fe55, 0.0);

    jac.set(P30, Fe56, 0.0);

    jac.set(P30, Co55, 0.0);

    jac.set(P30, Co56, 0.0);

    jac.set(P30, Co57, 0.0);

    jac.set(P30, Ni56, 0.0);

    jac.set(P30, Ni57, 0.0);

    jac.set(P30, Ni58, 0.0);

    jac.set(P30, Cu59, 0.0);

    jac.set(P30, Zn60, 0.0);

    jac.set(P31, N, 0.0);

    scratch = -screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(P31)*state.rho - screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho - screened_rates(k_p_P31_to_O16_O16_derived)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_He4_Al27_to_P31_reaclib)*Y(Al27)*state.rho - screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(P31)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    scratch = screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(Ne20)*state.rho;
    jac.set(P31, C12, scratch);

    jac.set(P31, C13, 0.0);

    jac.set(P31, C14, 0.0);

    jac.set(P31, N13, 0.0);

    jac.set(P31, N14, 0.0);

    jac.set(P31, N15, 0.0);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    jac.set(P31, O17, 0.0);

    jac.set(P31, O18, 0.0);

    jac.set(P31, F18, 0.0);

    jac.set(P31, F19, 0.0);

    scratch = screened_rates(k_C12_Ne20_to_p_P31_reaclib)*Y(C12)*state.rho;
    jac.set(P31, Ne20, scratch);

    jac.set(P31, Ne21, 0.0);

    jac.set(P31, Ne22, 0.0);

    jac.set(P31, Na22, 0.0);

    jac.set(P31, Na23, 0.0);

    jac.set(P31, Mg24, 0.0);

    jac.set(P31, Mg25, 0.0);

    jac.set(P31, Al26, 0.0);

    scratch = screened_rates(k_He4_Al27_to_P31_reaclib)*Y(He4)*state.rho;
    jac.set(P31, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    jac.set(P31, Si29, 0.0);

    jac.set(P31, P30, 0.0);

    scratch = -screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(He4)*state.rho - screened_rates(k_P31_to_He4_Al27_derived) - screened_rates(k_p_P31_to_C12_Ne20_derived)*Y(H1)*state.rho - screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho - screened_rates(k_p_P31_to_O16_O16_derived)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(P31, P31, scratch);

    scratch = screened_rates(k_S32_to_p_P31_derived);
    jac.set(P31, S32, scratch);

    jac.set(P31, S33, 0.0);

    jac.set(P31, Cl34, 0.0);

    scratch = screened_rates(k_Cl35_to_He4_P31_derived);
    jac.set(P31, Cl35, scratch);

    jac.set(P31, Ar36, 0.0);

    jac.set(P31, Ar37, 0.0);

    jac.set(P31, K38, 0.0);

    jac.set(P31, K39, 0.0);

    jac.set(P31, Ca40, 0.0);

    jac.set(P31, Ca41, 0.0);

    jac.set(P31, Sc42, 0.0);

    jac.set(P31, Sc43, 0.0);

    jac.set(P31, Ti44, 0.0);

    jac.set(P31, Ti45, 0.0);

    jac.set(P31, V46, 0.0);

    jac.set(P31, V47, 0.0);

    jac.set(P31, Cr48, 0.0);

    jac.set(P31, Cr49, 0.0);

    jac.set(P31, Mn50, 0.0);

    jac.set(P31, Mn51, 0.0);

    jac.set(P31, Mn52, 0.0);

    jac.set(P31, Fe52, 0.0);

    jac.set(P31, Fe53, 0.0);

    jac.set(P31, Fe54, 0.0);

    jac.set(P31, Fe55, 0.0);

    jac.set(P31, Fe56, 0.0);

    jac.set(P31, Co55, 0.0);

    jac.set(P31, Co56, 0.0);

    jac.set(P31, Co57, 0.0);

    jac.set(P31, Ni56, 0.0);

    jac.set(P31, Ni57, 0.0);

    jac.set(P31, Ni58, 0.0);

    jac.set(P31, Cu59, 0.0);

    jac.set(P31, Zn60, 0.0);

    jac.set(S32, N, 0.0);

    scratch = screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*state.rho + screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = -screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(S32)*state.rho - screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(S32)*state.rho + screened_rates(k_He4_Si28_to_S32_reaclib)*Y(Si28)*state.rho;
    jac.set(S32, He4, scratch);

    jac.set(S32, C12, 0.0);

    jac.set(S32, C13, 0.0);

    jac.set(S32, C14, 0.0);

    jac.set(S32, N13, 0.0);

    jac.set(S32, N14, 0.0);

    jac.set(S32, N15, 0.0);

    jac.set(S32, O16, 0.0);

    jac.set(S32, O17, 0.0);

    jac.set(S32, O18, 0.0);

    jac.set(S32, F18, 0.0);

    jac.set(S32, F19, 0.0);

    jac.set(S32, Ne20, 0.0);

    jac.set(S32, Ne21, 0.0);

    jac.set(S32, Ne22, 0.0);

    jac.set(S32, Na22, 0.0);

    jac.set(S32, Na23, 0.0);

    jac.set(S32, Mg24, 0.0);

    jac.set(S32, Mg25, 0.0);

    jac.set(S32, Al26, 0.0);

    jac.set(S32, Al27, 0.0);

    scratch = screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    jac.set(S32, Si29, 0.0);

    jac.set(S32, P30, 0.0);

    scratch = screened_rates(k_p_P31_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);

    scratch = -screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*state.rho - screened_rates(k_S32_to_He4_Si28_derived) - screened_rates(k_S32_to_p_P31_derived);
    jac.set(S32, S32, scratch);

    jac.set(S32, S33, 0.0);

    jac.set(S32, Cl34, 0.0);

    scratch = screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(H1)*state.rho;
    jac.set(S32, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_He4_S32_derived);
    jac.set(S32, Ar36, scratch);

    jac.set(S32, Ar37, 0.0);

    jac.set(S32, K38, 0.0);

    jac.set(S32, K39, 0.0);

    jac.set(S32, Ca40, 0.0);

    jac.set(S32, Ca41, 0.0);

    jac.set(S32, Sc42, 0.0);

    jac.set(S32, Sc43, 0.0);

    jac.set(S32, Ti44, 0.0);

    jac.set(S32, Ti45, 0.0);

    jac.set(S32, V46, 0.0);

    jac.set(S32, V47, 0.0);

    jac.set(S32, Cr48, 0.0);

    jac.set(S32, Cr49, 0.0);

    jac.set(S32, Mn50, 0.0);

    jac.set(S32, Mn51, 0.0);

    jac.set(S32, Mn52, 0.0);

    jac.set(S32, Fe52, 0.0);

    jac.set(S32, Fe53, 0.0);

    jac.set(S32, Fe54, 0.0);

    jac.set(S32, Fe55, 0.0);

    jac.set(S32, Fe56, 0.0);

    jac.set(S32, Co55, 0.0);

    jac.set(S32, Co56, 0.0);

    jac.set(S32, Co57, 0.0);

    jac.set(S32, Ni56, 0.0);

    jac.set(S32, Ni57, 0.0);

    jac.set(S32, Ni58, 0.0);

    jac.set(S32, Cu59, 0.0);

    jac.set(S32, Zn60, 0.0);

    jac.set(S33, N, 0.0);

    scratch = -screened_rates(k_p_S33_to_Cl34_reaclib)*Y(S33)*state.rho - screened_rates(k_p_S33_to_He4_P30_derived)*Y(S33)*state.rho;
    jac.set(S33, H1, scratch);

    scratch = screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(P30)*state.rho - screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(S33)*state.rho + screened_rates(k_He4_Si29_to_S33_reaclib)*Y(Si29)*state.rho;
    jac.set(S33, He4, scratch);

    jac.set(S33, C12, 0.0);

    jac.set(S33, C13, 0.0);

    jac.set(S33, C14, 0.0);

    jac.set(S33, N13, 0.0);

    jac.set(S33, N14, 0.0);

    jac.set(S33, N15, 0.0);

    jac.set(S33, O16, 0.0);

    jac.set(S33, O17, 0.0);

    jac.set(S33, O18, 0.0);

    jac.set(S33, F18, 0.0);

    jac.set(S33, F19, 0.0);

    jac.set(S33, Ne20, 0.0);

    jac.set(S33, Ne21, 0.0);

    jac.set(S33, Ne22, 0.0);

    jac.set(S33, Na22, 0.0);

    jac.set(S33, Na23, 0.0);

    jac.set(S33, Mg24, 0.0);

    jac.set(S33, Mg25, 0.0);

    jac.set(S33, Al26, 0.0);

    jac.set(S33, Al27, 0.0);

    jac.set(S33, Si28, 0.0);

    scratch = screened_rates(k_He4_Si29_to_S33_reaclib)*Y(He4)*state.rho;
    jac.set(S33, Si29, scratch);

    scratch = screened_rates(k_He4_P30_to_p_S33_reaclib)*Y(He4)*state.rho;
    jac.set(S33, P30, scratch);

    jac.set(S33, P31, 0.0);

    jac.set(S33, S32, 0.0);

    scratch = -screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(He4)*state.rho - screened_rates(k_S33_to_He4_Si29_derived) - screened_rates(k_p_S33_to_Cl34_reaclib)*Y(H1)*state.rho - screened_rates(k_p_S33_to_He4_P30_derived)*Y(H1)*state.rho;
    jac.set(S33, S33, scratch);

    scratch = screened_rates(k_Cl34_to_p_S33_derived);
    jac.set(S33, Cl34, scratch);

    jac.set(S33, Cl35, 0.0);

    jac.set(S33, Ar36, 0.0);

    scratch = screened_rates(k_Ar37_to_He4_S33_derived);
    jac.set(S33, Ar37, scratch);

    jac.set(S33, K38, 0.0);

    jac.set(S33, K39, 0.0);

    jac.set(S33, Ca40, 0.0);

    jac.set(S33, Ca41, 0.0);

    jac.set(S33, Sc42, 0.0);

    jac.set(S33, Sc43, 0.0);

    jac.set(S33, Ti44, 0.0);

    jac.set(S33, Ti45, 0.0);

    jac.set(S33, V46, 0.0);

    jac.set(S33, V47, 0.0);

    jac.set(S33, Cr48, 0.0);

    jac.set(S33, Cr49, 0.0);

    jac.set(S33, Mn50, 0.0);

    jac.set(S33, Mn51, 0.0);

    jac.set(S33, Mn52, 0.0);

    jac.set(S33, Fe52, 0.0);

    jac.set(S33, Fe53, 0.0);

    jac.set(S33, Fe54, 0.0);

    jac.set(S33, Fe55, 0.0);

    jac.set(S33, Fe56, 0.0);

    jac.set(S33, Co55, 0.0);

    jac.set(S33, Co56, 0.0);

    jac.set(S33, Co57, 0.0);

    jac.set(S33, Ni56, 0.0);

    jac.set(S33, Ni57, 0.0);

    jac.set(S33, Ni58, 0.0);

    jac.set(S33, Cu59, 0.0);

    jac.set(S33, Zn60, 0.0);

    jac.set(Cl34, N, 0.0);

    scratch = screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*state.rho + screened_rates(k_p_S33_to_Cl34_reaclib)*Y(S33)*state.rho;
    jac.set(Cl34, H1, scratch);

    scratch = -screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(Cl34)*state.rho - screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*state.rho + screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(P30)*state.rho;
    jac.set(Cl34, He4, scratch);

    jac.set(Cl34, C12, 0.0);

    jac.set(Cl34, C13, 0.0);

    jac.set(Cl34, C14, 0.0);

    jac.set(Cl34, N13, 0.0);

    jac.set(Cl34, N14, 0.0);

    jac.set(Cl34, N15, 0.0);

    jac.set(Cl34, O16, 0.0);

    jac.set(Cl34, O17, 0.0);

    jac.set(Cl34, O18, 0.0);

    jac.set(Cl34, F18, 0.0);

    jac.set(Cl34, F19, 0.0);

    jac.set(Cl34, Ne20, 0.0);

    jac.set(Cl34, Ne21, 0.0);

    jac.set(Cl34, Ne22, 0.0);

    jac.set(Cl34, Na22, 0.0);

    jac.set(Cl34, Na23, 0.0);

    jac.set(Cl34, Mg24, 0.0);

    jac.set(Cl34, Mg25, 0.0);

    jac.set(Cl34, Al26, 0.0);

    jac.set(Cl34, Al27, 0.0);

    jac.set(Cl34, Si28, 0.0);

    jac.set(Cl34, Si29, 0.0);

    scratch = screened_rates(k_He4_P30_to_Cl34_reaclib)*Y(He4)*state.rho;
    jac.set(Cl34, P30, scratch);

    jac.set(Cl34, P31, 0.0);

    jac.set(Cl34, S32, 0.0);

    scratch = screened_rates(k_p_S33_to_Cl34_reaclib)*Y(H1)*state.rho;
    jac.set(Cl34, S33, scratch);

    scratch = -screened_rates(k_Cl34_to_He4_P30_derived) - screened_rates(k_Cl34_to_p_S33_derived) - screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(He4)*state.rho;
    jac.set(Cl34, Cl34, scratch);

    jac.set(Cl34, Cl35, 0.0);

    jac.set(Cl34, Ar36, 0.0);

    scratch = screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(H1)*state.rho;
    jac.set(Cl34, Ar37, scratch);

    scratch = screened_rates(k_K38_to_He4_Cl34_derived);
    jac.set(Cl34, K38, scratch);

    jac.set(Cl34, K39, 0.0);

    jac.set(Cl34, Ca40, 0.0);

    jac.set(Cl34, Ca41, 0.0);

    jac.set(Cl34, Sc42, 0.0);

    jac.set(Cl34, Sc43, 0.0);

    jac.set(Cl34, Ti44, 0.0);

    jac.set(Cl34, Ti45, 0.0);

    jac.set(Cl34, V46, 0.0);

    jac.set(Cl34, V47, 0.0);

    jac.set(Cl34, Cr48, 0.0);

    jac.set(Cl34, Cr49, 0.0);

    jac.set(Cl34, Mn50, 0.0);

    jac.set(Cl34, Mn51, 0.0);

    jac.set(Cl34, Mn52, 0.0);

    jac.set(Cl34, Fe52, 0.0);

    jac.set(Cl34, Fe53, 0.0);

    jac.set(Cl34, Fe54, 0.0);

    jac.set(Cl34, Fe55, 0.0);

    jac.set(Cl34, Fe56, 0.0);

    jac.set(Cl34, Co55, 0.0);

    jac.set(Cl34, Co56, 0.0);

    jac.set(Cl34, Co57, 0.0);

    jac.set(Cl34, Ni56, 0.0);

    jac.set(Cl34, Ni57, 0.0);

    jac.set(Cl34, Ni58, 0.0);

    jac.set(Cl34, Cu59, 0.0);

    jac.set(Cl34, Zn60, 0.0);

    jac.set(Cl35, N, 0.0);

    scratch = -screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(Cl35)*state.rho - screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(Cl35)*state.rho;
    jac.set(Cl35, H1, scratch);

    scratch = -screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(Cl35)*state.rho + screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(P31)*state.rho + screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(S32)*state.rho;
    jac.set(Cl35, He4, scratch);

    jac.set(Cl35, C12, 0.0);

    jac.set(Cl35, C13, 0.0);

    jac.set(Cl35, C14, 0.0);

    jac.set(Cl35, N13, 0.0);

    jac.set(Cl35, N14, 0.0);

    jac.set(Cl35, N15, 0.0);

    jac.set(Cl35, O16, 0.0);

    jac.set(Cl35, O17, 0.0);

    jac.set(Cl35, O18, 0.0);

    jac.set(Cl35, F18, 0.0);

    jac.set(Cl35, F19, 0.0);

    jac.set(Cl35, Ne20, 0.0);

    jac.set(Cl35, Ne21, 0.0);

    jac.set(Cl35, Ne22, 0.0);

    jac.set(Cl35, Na22, 0.0);

    jac.set(Cl35, Na23, 0.0);

    jac.set(Cl35, Mg24, 0.0);

    jac.set(Cl35, Mg25, 0.0);

    jac.set(Cl35, Al26, 0.0);

    jac.set(Cl35, Al27, 0.0);

    jac.set(Cl35, Si28, 0.0);

    jac.set(Cl35, Si29, 0.0);

    jac.set(Cl35, P30, 0.0);

    scratch = screened_rates(k_He4_P31_to_Cl35_reaclib)*Y(He4)*state.rho;
    jac.set(Cl35, P31, scratch);

    scratch = screened_rates(k_He4_S32_to_p_Cl35_derived)*Y(He4)*state.rho;
    jac.set(Cl35, S32, scratch);

    jac.set(Cl35, S33, 0.0);

    jac.set(Cl35, Cl34, 0.0);

    scratch = -screened_rates(k_Cl35_to_He4_P31_derived) - screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Cl35_to_He4_S32_reaclib)*Y(H1)*state.rho;
    jac.set(Cl35, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_p_Cl35_derived);
    jac.set(Cl35, Ar36, scratch);

    jac.set(Cl35, Ar37, 0.0);

    jac.set(Cl35, K38, 0.0);

    scratch = screened_rates(k_K39_to_He4_Cl35_derived);
    jac.set(Cl35, K39, scratch);

    jac.set(Cl35, Ca40, 0.0);

    jac.set(Cl35, Ca41, 0.0);

    jac.set(Cl35, Sc42, 0.0);

    jac.set(Cl35, Sc43, 0.0);

    jac.set(Cl35, Ti44, 0.0);

    jac.set(Cl35, Ti45, 0.0);

    jac.set(Cl35, V46, 0.0);

    jac.set(Cl35, V47, 0.0);

    jac.set(Cl35, Cr48, 0.0);

    jac.set(Cl35, Cr49, 0.0);

    jac.set(Cl35, Mn50, 0.0);

    jac.set(Cl35, Mn51, 0.0);

    jac.set(Cl35, Mn52, 0.0);

    jac.set(Cl35, Fe52, 0.0);

    jac.set(Cl35, Fe53, 0.0);

    jac.set(Cl35, Fe54, 0.0);

    jac.set(Cl35, Fe55, 0.0);

    jac.set(Cl35, Fe56, 0.0);

    jac.set(Cl35, Co55, 0.0);

    jac.set(Cl35, Co56, 0.0);

    jac.set(Cl35, Co57, 0.0);

    jac.set(Cl35, Ni56, 0.0);

    jac.set(Cl35, Ni57, 0.0);

    jac.set(Cl35, Ni58, 0.0);

    jac.set(Cl35, Cu59, 0.0);

    jac.set(Cl35, Zn60, 0.0);

    jac.set(Ar36, N, 0.0);

    scratch = screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(Cl35)*state.rho + screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*state.rho;
    jac.set(Ar36, H1, scratch);

    scratch = -screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(Ar36)*state.rho - screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*state.rho + screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(S32)*state.rho;
    jac.set(Ar36, He4, scratch);

    jac.set(Ar36, C12, 0.0);

    jac.set(Ar36, C13, 0.0);

    jac.set(Ar36, C14, 0.0);

    jac.set(Ar36, N13, 0.0);

    jac.set(Ar36, N14, 0.0);

    jac.set(Ar36, N15, 0.0);

    jac.set(Ar36, O16, 0.0);

    jac.set(Ar36, O17, 0.0);

    jac.set(Ar36, O18, 0.0);

    jac.set(Ar36, F18, 0.0);

    jac.set(Ar36, F19, 0.0);

    jac.set(Ar36, Ne20, 0.0);

    jac.set(Ar36, Ne21, 0.0);

    jac.set(Ar36, Ne22, 0.0);

    jac.set(Ar36, Na22, 0.0);

    jac.set(Ar36, Na23, 0.0);

    jac.set(Ar36, Mg24, 0.0);

    jac.set(Ar36, Mg25, 0.0);

    jac.set(Ar36, Al26, 0.0);

    jac.set(Ar36, Al27, 0.0);

    jac.set(Ar36, Si28, 0.0);

    jac.set(Ar36, Si29, 0.0);

    jac.set(Ar36, P30, 0.0);

    jac.set(Ar36, P31, 0.0);

    scratch = screened_rates(k_He4_S32_to_Ar36_reaclib)*Y(He4)*state.rho;
    jac.set(Ar36, S32, scratch);

    jac.set(Ar36, S33, 0.0);

    jac.set(Ar36, Cl34, 0.0);

    scratch = screened_rates(k_p_Cl35_to_Ar36_reaclib)*Y(H1)*state.rho;
    jac.set(Ar36, Cl35, scratch);

    scratch = -screened_rates(k_Ar36_to_He4_S32_derived) - screened_rates(k_Ar36_to_p_Cl35_derived) - screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(He4)*state.rho;
    jac.set(Ar36, Ar36, scratch);

    jac.set(Ar36, Ar37, 0.0);

    jac.set(Ar36, K38, 0.0);

    scratch = screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(H1)*state.rho;
    jac.set(Ar36, K39, scratch);

    scratch = screened_rates(k_Ca40_to_He4_Ar36_derived);
    jac.set(Ar36, Ca40, scratch);

    jac.set(Ar36, Ca41, 0.0);

    jac.set(Ar36, Sc42, 0.0);

    jac.set(Ar36, Sc43, 0.0);

    jac.set(Ar36, Ti44, 0.0);

    jac.set(Ar36, Ti45, 0.0);

    jac.set(Ar36, V46, 0.0);

    jac.set(Ar36, V47, 0.0);

    jac.set(Ar36, Cr48, 0.0);

    jac.set(Ar36, Cr49, 0.0);

    jac.set(Ar36, Mn50, 0.0);

    jac.set(Ar36, Mn51, 0.0);

    jac.set(Ar36, Mn52, 0.0);

    jac.set(Ar36, Fe52, 0.0);

    jac.set(Ar36, Fe53, 0.0);

    jac.set(Ar36, Fe54, 0.0);

    jac.set(Ar36, Fe55, 0.0);

    jac.set(Ar36, Fe56, 0.0);

    jac.set(Ar36, Co55, 0.0);

    jac.set(Ar36, Co56, 0.0);

    jac.set(Ar36, Co57, 0.0);

    jac.set(Ar36, Ni56, 0.0);

    jac.set(Ar36, Ni57, 0.0);

    jac.set(Ar36, Ni58, 0.0);

    jac.set(Ar36, Cu59, 0.0);

    jac.set(Ar36, Zn60, 0.0);

    jac.set(Ar37, N, 0.0);

    scratch = -screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(Ar37)*state.rho - screened_rates(k_p_Ar37_to_K38_reaclib)*Y(Ar37)*state.rho;
    jac.set(Ar37, H1, scratch);

    scratch = -screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(Ar37)*state.rho + screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(Cl34)*state.rho + screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(S33)*state.rho;
    jac.set(Ar37, He4, scratch);

    jac.set(Ar37, C12, 0.0);

    jac.set(Ar37, C13, 0.0);

    jac.set(Ar37, C14, 0.0);

    jac.set(Ar37, N13, 0.0);

    jac.set(Ar37, N14, 0.0);

    jac.set(Ar37, N15, 0.0);

    jac.set(Ar37, O16, 0.0);

    jac.set(Ar37, O17, 0.0);

    jac.set(Ar37, O18, 0.0);

    jac.set(Ar37, F18, 0.0);

    jac.set(Ar37, F19, 0.0);

    jac.set(Ar37, Ne20, 0.0);

    jac.set(Ar37, Ne21, 0.0);

    jac.set(Ar37, Ne22, 0.0);

    jac.set(Ar37, Na22, 0.0);

    jac.set(Ar37, Na23, 0.0);

    jac.set(Ar37, Mg24, 0.0);

    jac.set(Ar37, Mg25, 0.0);

    jac.set(Ar37, Al26, 0.0);

    jac.set(Ar37, Al27, 0.0);

    jac.set(Ar37, Si28, 0.0);

    jac.set(Ar37, Si29, 0.0);

    jac.set(Ar37, P30, 0.0);

    jac.set(Ar37, P31, 0.0);

    jac.set(Ar37, S32, 0.0);

    scratch = screened_rates(k_He4_S33_to_Ar37_reaclib)*Y(He4)*state.rho;
    jac.set(Ar37, S33, scratch);

    scratch = screened_rates(k_He4_Cl34_to_p_Ar37_reaclib)*Y(He4)*state.rho;
    jac.set(Ar37, Cl34, scratch);

    jac.set(Ar37, Cl35, 0.0);

    jac.set(Ar37, Ar36, 0.0);

    scratch = -screened_rates(k_Ar37_to_He4_S33_derived) - screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Ar37_to_He4_Cl34_derived)*Y(H1)*state.rho - screened_rates(k_p_Ar37_to_K38_reaclib)*Y(H1)*state.rho;
    jac.set(Ar37, Ar37, scratch);

    scratch = screened_rates(k_K38_to_p_Ar37_derived);
    jac.set(Ar37, K38, scratch);

    jac.set(Ar37, K39, 0.0);

    jac.set(Ar37, Ca40, 0.0);

    scratch = screened_rates(k_Ca41_to_He4_Ar37_derived);
    jac.set(Ar37, Ca41, scratch);

    jac.set(Ar37, Sc42, 0.0);

    jac.set(Ar37, Sc43, 0.0);

    jac.set(Ar37, Ti44, 0.0);

    jac.set(Ar37, Ti45, 0.0);

    jac.set(Ar37, V46, 0.0);

    jac.set(Ar37, V47, 0.0);

    jac.set(Ar37, Cr48, 0.0);

    jac.set(Ar37, Cr49, 0.0);

    jac.set(Ar37, Mn50, 0.0);

    jac.set(Ar37, Mn51, 0.0);

    jac.set(Ar37, Mn52, 0.0);

    jac.set(Ar37, Fe52, 0.0);

    jac.set(Ar37, Fe53, 0.0);

    jac.set(Ar37, Fe54, 0.0);

    jac.set(Ar37, Fe55, 0.0);

    jac.set(Ar37, Fe56, 0.0);

    jac.set(Ar37, Co55, 0.0);

    jac.set(Ar37, Co56, 0.0);

    jac.set(Ar37, Co57, 0.0);

    jac.set(Ar37, Ni56, 0.0);

    jac.set(Ar37, Ni57, 0.0);

    jac.set(Ar37, Ni58, 0.0);

    jac.set(Ar37, Cu59, 0.0);

    jac.set(Ar37, Zn60, 0.0);

    jac.set(K38, N, 0.0);

    scratch = screened_rates(k_p_Ar37_to_K38_reaclib)*Y(Ar37)*state.rho + screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*state.rho;
    jac.set(K38, H1, scratch);

    scratch = screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(Cl34)*state.rho - screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(K38)*state.rho - screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(K38)*state.rho;
    jac.set(K38, He4, scratch);

    jac.set(K38, C12, 0.0);

    jac.set(K38, C13, 0.0);

    jac.set(K38, C14, 0.0);

    jac.set(K38, N13, 0.0);

    jac.set(K38, N14, 0.0);

    jac.set(K38, N15, 0.0);

    jac.set(K38, O16, 0.0);

    jac.set(K38, O17, 0.0);

    jac.set(K38, O18, 0.0);

    jac.set(K38, F18, 0.0);

    jac.set(K38, F19, 0.0);

    jac.set(K38, Ne20, 0.0);

    jac.set(K38, Ne21, 0.0);

    jac.set(K38, Ne22, 0.0);

    jac.set(K38, Na22, 0.0);

    jac.set(K38, Na23, 0.0);

    jac.set(K38, Mg24, 0.0);

    jac.set(K38, Mg25, 0.0);

    jac.set(K38, Al26, 0.0);

    jac.set(K38, Al27, 0.0);

    jac.set(K38, Si28, 0.0);

    jac.set(K38, Si29, 0.0);

    jac.set(K38, P30, 0.0);

    jac.set(K38, P31, 0.0);

    jac.set(K38, S32, 0.0);

    jac.set(K38, S33, 0.0);

    scratch = screened_rates(k_He4_Cl34_to_K38_reaclib)*Y(He4)*state.rho;
    jac.set(K38, Cl34, scratch);

    jac.set(K38, Cl35, 0.0);

    jac.set(K38, Ar36, 0.0);

    scratch = screened_rates(k_p_Ar37_to_K38_reaclib)*Y(H1)*state.rho;
    jac.set(K38, Ar37, scratch);

    scratch = -screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*state.rho - screened_rates(k_K38_to_He4_Cl34_derived) - screened_rates(k_K38_to_p_Ar37_derived);
    jac.set(K38, K38, scratch);

    jac.set(K38, K39, 0.0);

    jac.set(K38, Ca40, 0.0);

    scratch = screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(H1)*state.rho;
    jac.set(K38, Ca41, scratch);

    scratch = screened_rates(k_Sc42_to_He4_K38_derived);
    jac.set(K38, Sc42, scratch);

    jac.set(K38, Sc43, 0.0);

    jac.set(K38, Ti44, 0.0);

    jac.set(K38, Ti45, 0.0);

    jac.set(K38, V46, 0.0);

    jac.set(K38, V47, 0.0);

    jac.set(K38, Cr48, 0.0);

    jac.set(K38, Cr49, 0.0);

    jac.set(K38, Mn50, 0.0);

    jac.set(K38, Mn51, 0.0);

    jac.set(K38, Mn52, 0.0);

    jac.set(K38, Fe52, 0.0);

    jac.set(K38, Fe53, 0.0);

    jac.set(K38, Fe54, 0.0);

    jac.set(K38, Fe55, 0.0);

    jac.set(K38, Fe56, 0.0);

    jac.set(K38, Co55, 0.0);

    jac.set(K38, Co56, 0.0);

    jac.set(K38, Co57, 0.0);

    jac.set(K38, Ni56, 0.0);

    jac.set(K38, Ni57, 0.0);

    jac.set(K38, Ni58, 0.0);

    jac.set(K38, Cu59, 0.0);

    jac.set(K38, Zn60, 0.0);

    jac.set(K39, N, 0.0);

    scratch = -screened_rates(k_p_K39_to_Ca40_reaclib)*Y(K39)*state.rho - screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(K39)*state.rho;
    jac.set(K39, H1, scratch);

    scratch = screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(Ar36)*state.rho + screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(Cl35)*state.rho - screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(K39)*state.rho;
    jac.set(K39, He4, scratch);

    jac.set(K39, C12, 0.0);

    jac.set(K39, C13, 0.0);

    jac.set(K39, C14, 0.0);

    jac.set(K39, N13, 0.0);

    jac.set(K39, N14, 0.0);

    jac.set(K39, N15, 0.0);

    jac.set(K39, O16, 0.0);

    jac.set(K39, O17, 0.0);

    jac.set(K39, O18, 0.0);

    jac.set(K39, F18, 0.0);

    jac.set(K39, F19, 0.0);

    jac.set(K39, Ne20, 0.0);

    jac.set(K39, Ne21, 0.0);

    jac.set(K39, Ne22, 0.0);

    jac.set(K39, Na22, 0.0);

    jac.set(K39, Na23, 0.0);

    jac.set(K39, Mg24, 0.0);

    jac.set(K39, Mg25, 0.0);

    jac.set(K39, Al26, 0.0);

    jac.set(K39, Al27, 0.0);

    jac.set(K39, Si28, 0.0);

    jac.set(K39, Si29, 0.0);

    jac.set(K39, P30, 0.0);

    jac.set(K39, P31, 0.0);

    jac.set(K39, S32, 0.0);

    jac.set(K39, S33, 0.0);

    jac.set(K39, Cl34, 0.0);

    scratch = screened_rates(k_He4_Cl35_to_K39_reaclib)*Y(He4)*state.rho;
    jac.set(K39, Cl35, scratch);

    scratch = screened_rates(k_He4_Ar36_to_p_K39_derived)*Y(He4)*state.rho;
    jac.set(K39, Ar36, scratch);

    jac.set(K39, Ar37, 0.0);

    jac.set(K39, K38, 0.0);

    scratch = -screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(He4)*state.rho - screened_rates(k_K39_to_He4_Cl35_derived) - screened_rates(k_p_K39_to_Ca40_reaclib)*Y(H1)*state.rho - screened_rates(k_p_K39_to_He4_Ar36_reaclib)*Y(H1)*state.rho;
    jac.set(K39, K39, scratch);

    scratch = screened_rates(k_Ca40_to_p_K39_derived);
    jac.set(K39, Ca40, scratch);

    jac.set(K39, Ca41, 0.0);

    jac.set(K39, Sc42, 0.0);

    scratch = screened_rates(k_Sc43_to_He4_K39_derived);
    jac.set(K39, Sc43, scratch);

    jac.set(K39, Ti44, 0.0);

    jac.set(K39, Ti45, 0.0);

    jac.set(K39, V46, 0.0);

    jac.set(K39, V47, 0.0);

    jac.set(K39, Cr48, 0.0);

    jac.set(K39, Cr49, 0.0);

    jac.set(K39, Mn50, 0.0);

    jac.set(K39, Mn51, 0.0);

    jac.set(K39, Mn52, 0.0);

    jac.set(K39, Fe52, 0.0);

    jac.set(K39, Fe53, 0.0);

    jac.set(K39, Fe54, 0.0);

    jac.set(K39, Fe55, 0.0);

    jac.set(K39, Fe56, 0.0);

    jac.set(K39, Co55, 0.0);

    jac.set(K39, Co56, 0.0);

    jac.set(K39, Co57, 0.0);

    jac.set(K39, Ni56, 0.0);

    jac.set(K39, Ni57, 0.0);

    jac.set(K39, Ni58, 0.0);

    jac.set(K39, Cu59, 0.0);

    jac.set(K39, Zn60, 0.0);

    jac.set(Ca40, N, 0.0);

    scratch = screened_rates(k_p_K39_to_Ca40_reaclib)*Y(K39)*state.rho + screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho;
    jac.set(Ca40, H1, scratch);

    scratch = screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(Ar36)*state.rho - screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(Ca40)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*state.rho;
    jac.set(Ca40, He4, scratch);

    jac.set(Ca40, C12, 0.0);

    jac.set(Ca40, C13, 0.0);

    jac.set(Ca40, C14, 0.0);

    jac.set(Ca40, N13, 0.0);

    jac.set(Ca40, N14, 0.0);

    jac.set(Ca40, N15, 0.0);

    jac.set(Ca40, O16, 0.0);

    jac.set(Ca40, O17, 0.0);

    jac.set(Ca40, O18, 0.0);

    jac.set(Ca40, F18, 0.0);

    jac.set(Ca40, F19, 0.0);

    jac.set(Ca40, Ne20, 0.0);

    jac.set(Ca40, Ne21, 0.0);

    jac.set(Ca40, Ne22, 0.0);

    jac.set(Ca40, Na22, 0.0);

    jac.set(Ca40, Na23, 0.0);

    jac.set(Ca40, Mg24, 0.0);

    jac.set(Ca40, Mg25, 0.0);

    jac.set(Ca40, Al26, 0.0);

    jac.set(Ca40, Al27, 0.0);

    jac.set(Ca40, Si28, 0.0);

    jac.set(Ca40, Si29, 0.0);

    jac.set(Ca40, P30, 0.0);

    jac.set(Ca40, P31, 0.0);

    jac.set(Ca40, S32, 0.0);

    jac.set(Ca40, S33, 0.0);

    jac.set(Ca40, Cl34, 0.0);

    jac.set(Ca40, Cl35, 0.0);

    scratch = screened_rates(k_He4_Ar36_to_Ca40_reaclib)*Y(He4)*state.rho;
    jac.set(Ca40, Ar36, scratch);

    jac.set(Ca40, Ar37, 0.0);

    jac.set(Ca40, K38, 0.0);

    scratch = screened_rates(k_p_K39_to_Ca40_reaclib)*Y(H1)*state.rho;
    jac.set(Ca40, K39, scratch);

    scratch = -screened_rates(k_Ca40_to_He4_Ar36_derived) - screened_rates(k_Ca40_to_p_K39_derived) - screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(He4)*state.rho;
    jac.set(Ca40, Ca40, scratch);

    jac.set(Ca40, Ca41, 0.0);

    jac.set(Ca40, Sc42, 0.0);

    scratch = screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(H1)*state.rho;
    jac.set(Ca40, Sc43, scratch);

    scratch = screened_rates(k_Ti44_to_He4_Ca40_derived);
    jac.set(Ca40, Ti44, scratch);

    jac.set(Ca40, Ti45, 0.0);

    jac.set(Ca40, V46, 0.0);

    jac.set(Ca40, V47, 0.0);

    jac.set(Ca40, Cr48, 0.0);

    jac.set(Ca40, Cr49, 0.0);

    jac.set(Ca40, Mn50, 0.0);

    jac.set(Ca40, Mn51, 0.0);

    jac.set(Ca40, Mn52, 0.0);

    jac.set(Ca40, Fe52, 0.0);

    jac.set(Ca40, Fe53, 0.0);

    jac.set(Ca40, Fe54, 0.0);

    jac.set(Ca40, Fe55, 0.0);

    jac.set(Ca40, Fe56, 0.0);

    jac.set(Ca40, Co55, 0.0);

    jac.set(Ca40, Co56, 0.0);

    jac.set(Ca40, Co57, 0.0);

    jac.set(Ca40, Ni56, 0.0);

    jac.set(Ca40, Ni57, 0.0);

    jac.set(Ca40, Ni58, 0.0);

    jac.set(Ca40, Cu59, 0.0);

    jac.set(Ca40, Zn60, 0.0);

    jac.set(Ca41, N, 0.0);

    scratch = -screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(Ca41)*state.rho - screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(Ca41)*state.rho;
    jac.set(Ca41, H1, scratch);

    scratch = screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(Ar37)*state.rho - screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(Ca41)*state.rho + screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(K38)*state.rho;
    jac.set(Ca41, He4, scratch);

    jac.set(Ca41, C12, 0.0);

    jac.set(Ca41, C13, 0.0);

    jac.set(Ca41, C14, 0.0);

    jac.set(Ca41, N13, 0.0);

    jac.set(Ca41, N14, 0.0);

    jac.set(Ca41, N15, 0.0);

    jac.set(Ca41, O16, 0.0);

    jac.set(Ca41, O17, 0.0);

    jac.set(Ca41, O18, 0.0);

    jac.set(Ca41, F18, 0.0);

    jac.set(Ca41, F19, 0.0);

    jac.set(Ca41, Ne20, 0.0);

    jac.set(Ca41, Ne21, 0.0);

    jac.set(Ca41, Ne22, 0.0);

    jac.set(Ca41, Na22, 0.0);

    jac.set(Ca41, Na23, 0.0);

    jac.set(Ca41, Mg24, 0.0);

    jac.set(Ca41, Mg25, 0.0);

    jac.set(Ca41, Al26, 0.0);

    jac.set(Ca41, Al27, 0.0);

    jac.set(Ca41, Si28, 0.0);

    jac.set(Ca41, Si29, 0.0);

    jac.set(Ca41, P30, 0.0);

    jac.set(Ca41, P31, 0.0);

    jac.set(Ca41, S32, 0.0);

    jac.set(Ca41, S33, 0.0);

    jac.set(Ca41, Cl34, 0.0);

    jac.set(Ca41, Cl35, 0.0);

    jac.set(Ca41, Ar36, 0.0);

    scratch = screened_rates(k_He4_Ar37_to_Ca41_reaclib)*Y(He4)*state.rho;
    jac.set(Ca41, Ar37, scratch);

    scratch = screened_rates(k_He4_K38_to_p_Ca41_reaclib)*Y(He4)*state.rho;
    jac.set(Ca41, K38, scratch);

    jac.set(Ca41, K39, 0.0);

    jac.set(Ca41, Ca40, 0.0);

    scratch = -screened_rates(k_Ca41_to_He4_Ar37_derived) - screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Ca41_to_He4_K38_derived)*Y(H1)*state.rho - screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(H1)*state.rho;
    jac.set(Ca41, Ca41, scratch);

    scratch = screened_rates(k_Sc42_to_p_Ca41_derived);
    jac.set(Ca41, Sc42, scratch);

    jac.set(Ca41, Sc43, 0.0);

    jac.set(Ca41, Ti44, 0.0);

    scratch = screened_rates(k_Ti45_to_He4_Ca41_derived);
    jac.set(Ca41, Ti45, scratch);

    jac.set(Ca41, V46, 0.0);

    jac.set(Ca41, V47, 0.0);

    jac.set(Ca41, Cr48, 0.0);

    jac.set(Ca41, Cr49, 0.0);

    jac.set(Ca41, Mn50, 0.0);

    jac.set(Ca41, Mn51, 0.0);

    jac.set(Ca41, Mn52, 0.0);

    jac.set(Ca41, Fe52, 0.0);

    jac.set(Ca41, Fe53, 0.0);

    jac.set(Ca41, Fe54, 0.0);

    jac.set(Ca41, Fe55, 0.0);

    jac.set(Ca41, Fe56, 0.0);

    jac.set(Ca41, Co55, 0.0);

    jac.set(Ca41, Co56, 0.0);

    jac.set(Ca41, Co57, 0.0);

    jac.set(Ca41, Ni56, 0.0);

    jac.set(Ca41, Ni57, 0.0);

    jac.set(Ca41, Ni58, 0.0);

    jac.set(Ca41, Cu59, 0.0);

    jac.set(Ca41, Zn60, 0.0);

    jac.set(Sc42, N, 0.0);

    scratch = screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(Ca41)*state.rho + screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*state.rho;
    jac.set(Sc42, H1, scratch);

    scratch = screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(K38)*state.rho - screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(Sc42)*state.rho - screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(Sc42)*state.rho;
    jac.set(Sc42, He4, scratch);

    jac.set(Sc42, C12, 0.0);

    jac.set(Sc42, C13, 0.0);

    jac.set(Sc42, C14, 0.0);

    jac.set(Sc42, N13, 0.0);

    jac.set(Sc42, N14, 0.0);

    jac.set(Sc42, N15, 0.0);

    jac.set(Sc42, O16, 0.0);

    jac.set(Sc42, O17, 0.0);

    jac.set(Sc42, O18, 0.0);

    jac.set(Sc42, F18, 0.0);

    jac.set(Sc42, F19, 0.0);

    jac.set(Sc42, Ne20, 0.0);

    jac.set(Sc42, Ne21, 0.0);

    jac.set(Sc42, Ne22, 0.0);

    jac.set(Sc42, Na22, 0.0);

    jac.set(Sc42, Na23, 0.0);

    jac.set(Sc42, Mg24, 0.0);

    jac.set(Sc42, Mg25, 0.0);

    jac.set(Sc42, Al26, 0.0);

    jac.set(Sc42, Al27, 0.0);

    jac.set(Sc42, Si28, 0.0);

    jac.set(Sc42, Si29, 0.0);

    jac.set(Sc42, P30, 0.0);

    jac.set(Sc42, P31, 0.0);

    jac.set(Sc42, S32, 0.0);

    jac.set(Sc42, S33, 0.0);

    jac.set(Sc42, Cl34, 0.0);

    jac.set(Sc42, Cl35, 0.0);

    jac.set(Sc42, Ar36, 0.0);

    jac.set(Sc42, Ar37, 0.0);

    scratch = screened_rates(k_He4_K38_to_Sc42_reaclib)*Y(He4)*state.rho;
    jac.set(Sc42, K38, scratch);

    jac.set(Sc42, K39, 0.0);

    jac.set(Sc42, Ca40, 0.0);

    scratch = screened_rates(k_p_Ca41_to_Sc42_reaclib)*Y(H1)*state.rho;
    jac.set(Sc42, Ca41, scratch);

    scratch = -screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*state.rho - screened_rates(k_Sc42_to_He4_K38_derived) - screened_rates(k_Sc42_to_p_Ca41_derived);
    jac.set(Sc42, Sc42, scratch);

    jac.set(Sc42, Sc43, 0.0);

    jac.set(Sc42, Ti44, 0.0);

    scratch = screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(H1)*state.rho;
    jac.set(Sc42, Ti45, scratch);

    scratch = screened_rates(k_V46_to_He4_Sc42_derived);
    jac.set(Sc42, V46, scratch);

    jac.set(Sc42, V47, 0.0);

    jac.set(Sc42, Cr48, 0.0);

    jac.set(Sc42, Cr49, 0.0);

    jac.set(Sc42, Mn50, 0.0);

    jac.set(Sc42, Mn51, 0.0);

    jac.set(Sc42, Mn52, 0.0);

    jac.set(Sc42, Fe52, 0.0);

    jac.set(Sc42, Fe53, 0.0);

    jac.set(Sc42, Fe54, 0.0);

    jac.set(Sc42, Fe55, 0.0);

    jac.set(Sc42, Fe56, 0.0);

    jac.set(Sc42, Co55, 0.0);

    jac.set(Sc42, Co56, 0.0);

    jac.set(Sc42, Co57, 0.0);

    jac.set(Sc42, Ni56, 0.0);

    jac.set(Sc42, Ni57, 0.0);

    jac.set(Sc42, Ni58, 0.0);

    jac.set(Sc42, Cu59, 0.0);

    jac.set(Sc42, Zn60, 0.0);

    jac.set(Sc43, N, 0.0);

    scratch = -screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho - screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(Sc43)*state.rho;
    jac.set(Sc43, H1, scratch);

    scratch = screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(Ca40)*state.rho + screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(K39)*state.rho - screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(Sc43)*state.rho;
    jac.set(Sc43, He4, scratch);

    jac.set(Sc43, C12, 0.0);

    jac.set(Sc43, C13, 0.0);

    jac.set(Sc43, C14, 0.0);

    jac.set(Sc43, N13, 0.0);

    jac.set(Sc43, N14, 0.0);

    jac.set(Sc43, N15, 0.0);

    jac.set(Sc43, O16, 0.0);

    jac.set(Sc43, O17, 0.0);

    jac.set(Sc43, O18, 0.0);

    jac.set(Sc43, F18, 0.0);

    jac.set(Sc43, F19, 0.0);

    jac.set(Sc43, Ne20, 0.0);

    jac.set(Sc43, Ne21, 0.0);

    jac.set(Sc43, Ne22, 0.0);

    jac.set(Sc43, Na22, 0.0);

    jac.set(Sc43, Na23, 0.0);

    jac.set(Sc43, Mg24, 0.0);

    jac.set(Sc43, Mg25, 0.0);

    jac.set(Sc43, Al26, 0.0);

    jac.set(Sc43, Al27, 0.0);

    jac.set(Sc43, Si28, 0.0);

    jac.set(Sc43, Si29, 0.0);

    jac.set(Sc43, P30, 0.0);

    jac.set(Sc43, P31, 0.0);

    jac.set(Sc43, S32, 0.0);

    jac.set(Sc43, S33, 0.0);

    jac.set(Sc43, Cl34, 0.0);

    jac.set(Sc43, Cl35, 0.0);

    jac.set(Sc43, Ar36, 0.0);

    jac.set(Sc43, Ar37, 0.0);

    jac.set(Sc43, K38, 0.0);

    scratch = screened_rates(k_He4_K39_to_Sc43_reaclib)*Y(He4)*state.rho;
    jac.set(Sc43, K39, scratch);

    scratch = screened_rates(k_He4_Ca40_to_p_Sc43_derived)*Y(He4)*state.rho;
    jac.set(Sc43, Ca40, scratch);

    jac.set(Sc43, Ca41, 0.0);

    jac.set(Sc43, Sc42, 0.0);

    scratch = -screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(He4)*state.rho - screened_rates(k_Sc43_to_He4_K39_derived) - screened_rates(k_p_Sc43_to_He4_Ca40_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(H1)*state.rho;
    jac.set(Sc43, Sc43, scratch);

    scratch = screened_rates(k_Ti44_to_p_Sc43_derived);
    jac.set(Sc43, Ti44, scratch);

    jac.set(Sc43, Ti45, 0.0);

    jac.set(Sc43, V46, 0.0);

    scratch = screened_rates(k_V47_to_He4_Sc43_derived);
    jac.set(Sc43, V47, scratch);

    jac.set(Sc43, Cr48, 0.0);

    jac.set(Sc43, Cr49, 0.0);

    jac.set(Sc43, Mn50, 0.0);

    jac.set(Sc43, Mn51, 0.0);

    jac.set(Sc43, Mn52, 0.0);

    jac.set(Sc43, Fe52, 0.0);

    jac.set(Sc43, Fe53, 0.0);

    jac.set(Sc43, Fe54, 0.0);

    jac.set(Sc43, Fe55, 0.0);

    jac.set(Sc43, Fe56, 0.0);

    jac.set(Sc43, Co55, 0.0);

    jac.set(Sc43, Co56, 0.0);

    jac.set(Sc43, Co57, 0.0);

    jac.set(Sc43, Ni56, 0.0);

    jac.set(Sc43, Ni57, 0.0);

    jac.set(Sc43, Ni58, 0.0);

    jac.set(Sc43, Cu59, 0.0);

    jac.set(Sc43, Zn60, 0.0);

    jac.set(Ti44, N, 0.0);

    scratch = screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(Sc43)*state.rho + screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*state.rho;
    jac.set(Ti44, H1, scratch);

    scratch = screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(Ca40)*state.rho - screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(Ti44)*state.rho - screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(Ti44)*state.rho;
    jac.set(Ti44, He4, scratch);

    jac.set(Ti44, C12, 0.0);

    jac.set(Ti44, C13, 0.0);

    jac.set(Ti44, C14, 0.0);

    jac.set(Ti44, N13, 0.0);

    jac.set(Ti44, N14, 0.0);

    jac.set(Ti44, N15, 0.0);

    jac.set(Ti44, O16, 0.0);

    jac.set(Ti44, O17, 0.0);

    jac.set(Ti44, O18, 0.0);

    jac.set(Ti44, F18, 0.0);

    jac.set(Ti44, F19, 0.0);

    jac.set(Ti44, Ne20, 0.0);

    jac.set(Ti44, Ne21, 0.0);

    jac.set(Ti44, Ne22, 0.0);

    jac.set(Ti44, Na22, 0.0);

    jac.set(Ti44, Na23, 0.0);

    jac.set(Ti44, Mg24, 0.0);

    jac.set(Ti44, Mg25, 0.0);

    jac.set(Ti44, Al26, 0.0);

    jac.set(Ti44, Al27, 0.0);

    jac.set(Ti44, Si28, 0.0);

    jac.set(Ti44, Si29, 0.0);

    jac.set(Ti44, P30, 0.0);

    jac.set(Ti44, P31, 0.0);

    jac.set(Ti44, S32, 0.0);

    jac.set(Ti44, S33, 0.0);

    jac.set(Ti44, Cl34, 0.0);

    jac.set(Ti44, Cl35, 0.0);

    jac.set(Ti44, Ar36, 0.0);

    jac.set(Ti44, Ar37, 0.0);

    jac.set(Ti44, K38, 0.0);

    jac.set(Ti44, K39, 0.0);

    scratch = screened_rates(k_He4_Ca40_to_Ti44_reaclib)*Y(He4)*state.rho;
    jac.set(Ti44, Ca40, scratch);

    jac.set(Ti44, Ca41, 0.0);

    jac.set(Ti44, Sc42, 0.0);

    scratch = screened_rates(k_p_Sc43_to_Ti44_reaclib)*Y(H1)*state.rho;
    jac.set(Ti44, Sc43, scratch);

    scratch = -screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*state.rho - screened_rates(k_Ti44_to_He4_Ca40_derived) - screened_rates(k_Ti44_to_p_Sc43_derived);
    jac.set(Ti44, Ti44, scratch);

    jac.set(Ti44, Ti45, 0.0);

    jac.set(Ti44, V46, 0.0);

    scratch = screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(H1)*state.rho;
    jac.set(Ti44, V47, scratch);

    scratch = screened_rates(k_Cr48_to_He4_Ti44_derived);
    jac.set(Ti44, Cr48, scratch);

    jac.set(Ti44, Cr49, 0.0);

    jac.set(Ti44, Mn50, 0.0);

    jac.set(Ti44, Mn51, 0.0);

    jac.set(Ti44, Mn52, 0.0);

    jac.set(Ti44, Fe52, 0.0);

    jac.set(Ti44, Fe53, 0.0);

    jac.set(Ti44, Fe54, 0.0);

    jac.set(Ti44, Fe55, 0.0);

    jac.set(Ti44, Fe56, 0.0);

    jac.set(Ti44, Co55, 0.0);

    jac.set(Ti44, Co56, 0.0);

    jac.set(Ti44, Co57, 0.0);

    jac.set(Ti44, Ni56, 0.0);

    jac.set(Ti44, Ni57, 0.0);

    jac.set(Ti44, Ni58, 0.0);

    jac.set(Ti44, Cu59, 0.0);

    jac.set(Ti44, Zn60, 0.0);

    jac.set(Ti45, N, 0.0);

    scratch = -screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(Ti45)*state.rho - screened_rates(k_p_Ti45_to_V46_reaclib)*Y(Ti45)*state.rho;
    jac.set(Ti45, H1, scratch);

    scratch = screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(Ca41)*state.rho + screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(Sc42)*state.rho - screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(Ti45)*state.rho;
    jac.set(Ti45, He4, scratch);

    jac.set(Ti45, C12, 0.0);

    jac.set(Ti45, C13, 0.0);

    jac.set(Ti45, C14, 0.0);

    jac.set(Ti45, N13, 0.0);

    jac.set(Ti45, N14, 0.0);

    jac.set(Ti45, N15, 0.0);

    jac.set(Ti45, O16, 0.0);

    jac.set(Ti45, O17, 0.0);

    jac.set(Ti45, O18, 0.0);

    jac.set(Ti45, F18, 0.0);

    jac.set(Ti45, F19, 0.0);

    jac.set(Ti45, Ne20, 0.0);

    jac.set(Ti45, Ne21, 0.0);

    jac.set(Ti45, Ne22, 0.0);

    jac.set(Ti45, Na22, 0.0);

    jac.set(Ti45, Na23, 0.0);

    jac.set(Ti45, Mg24, 0.0);

    jac.set(Ti45, Mg25, 0.0);

    jac.set(Ti45, Al26, 0.0);

    jac.set(Ti45, Al27, 0.0);

    jac.set(Ti45, Si28, 0.0);

    jac.set(Ti45, Si29, 0.0);

    jac.set(Ti45, P30, 0.0);

    jac.set(Ti45, P31, 0.0);

    jac.set(Ti45, S32, 0.0);

    jac.set(Ti45, S33, 0.0);

    jac.set(Ti45, Cl34, 0.0);

    jac.set(Ti45, Cl35, 0.0);

    jac.set(Ti45, Ar36, 0.0);

    jac.set(Ti45, Ar37, 0.0);

    jac.set(Ti45, K38, 0.0);

    jac.set(Ti45, K39, 0.0);

    jac.set(Ti45, Ca40, 0.0);

    scratch = screened_rates(k_He4_Ca41_to_Ti45_reaclib)*Y(He4)*state.rho;
    jac.set(Ti45, Ca41, scratch);

    scratch = screened_rates(k_He4_Sc42_to_p_Ti45_reaclib)*Y(He4)*state.rho;
    jac.set(Ti45, Sc42, scratch);

    jac.set(Ti45, Sc43, 0.0);

    jac.set(Ti45, Ti44, 0.0);

    scratch = -screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(He4)*state.rho - screened_rates(k_Ti45_to_He4_Ca41_derived) - screened_rates(k_p_Ti45_to_He4_Sc42_derived)*Y(H1)*state.rho - screened_rates(k_p_Ti45_to_V46_reaclib)*Y(H1)*state.rho;
    jac.set(Ti45, Ti45, scratch);

    scratch = screened_rates(k_V46_to_p_Ti45_derived);
    jac.set(Ti45, V46, scratch);

    jac.set(Ti45, V47, 0.0);

    jac.set(Ti45, Cr48, 0.0);

    scratch = screened_rates(k_Cr49_to_He4_Ti45_derived);
    jac.set(Ti45, Cr49, scratch);

    jac.set(Ti45, Mn50, 0.0);

    jac.set(Ti45, Mn51, 0.0);

    jac.set(Ti45, Mn52, 0.0);

    jac.set(Ti45, Fe52, 0.0);

    jac.set(Ti45, Fe53, 0.0);

    jac.set(Ti45, Fe54, 0.0);

    jac.set(Ti45, Fe55, 0.0);

    jac.set(Ti45, Fe56, 0.0);

    jac.set(Ti45, Co55, 0.0);

    jac.set(Ti45, Co56, 0.0);

    jac.set(Ti45, Co57, 0.0);

    jac.set(Ti45, Ni56, 0.0);

    jac.set(Ti45, Ni57, 0.0);

    jac.set(Ti45, Ni58, 0.0);

    jac.set(Ti45, Cu59, 0.0);

    jac.set(Ti45, Zn60, 0.0);

    jac.set(V46, N, 0.0);

    scratch = screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*state.rho + screened_rates(k_p_Ti45_to_V46_reaclib)*Y(Ti45)*state.rho;
    jac.set(V46, H1, scratch);

    scratch = screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(Sc42)*state.rho - screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(V46)*state.rho - screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(V46)*state.rho;
    jac.set(V46, He4, scratch);

    jac.set(V46, C12, 0.0);

    jac.set(V46, C13, 0.0);

    jac.set(V46, C14, 0.0);

    jac.set(V46, N13, 0.0);

    jac.set(V46, N14, 0.0);

    jac.set(V46, N15, 0.0);

    jac.set(V46, O16, 0.0);

    jac.set(V46, O17, 0.0);

    jac.set(V46, O18, 0.0);

    jac.set(V46, F18, 0.0);

    jac.set(V46, F19, 0.0);

    jac.set(V46, Ne20, 0.0);

    jac.set(V46, Ne21, 0.0);

    jac.set(V46, Ne22, 0.0);

    jac.set(V46, Na22, 0.0);

    jac.set(V46, Na23, 0.0);

    jac.set(V46, Mg24, 0.0);

    jac.set(V46, Mg25, 0.0);

    jac.set(V46, Al26, 0.0);

    jac.set(V46, Al27, 0.0);

    jac.set(V46, Si28, 0.0);

    jac.set(V46, Si29, 0.0);

    jac.set(V46, P30, 0.0);

    jac.set(V46, P31, 0.0);

    jac.set(V46, S32, 0.0);

    jac.set(V46, S33, 0.0);

    jac.set(V46, Cl34, 0.0);

    jac.set(V46, Cl35, 0.0);

    jac.set(V46, Ar36, 0.0);

    jac.set(V46, Ar37, 0.0);

    jac.set(V46, K38, 0.0);

    jac.set(V46, K39, 0.0);

    jac.set(V46, Ca40, 0.0);

    jac.set(V46, Ca41, 0.0);

    scratch = screened_rates(k_He4_Sc42_to_V46_reaclib)*Y(He4)*state.rho;
    jac.set(V46, Sc42, scratch);

    jac.set(V46, Sc43, 0.0);

    jac.set(V46, Ti44, 0.0);

    scratch = screened_rates(k_p_Ti45_to_V46_reaclib)*Y(H1)*state.rho;
    jac.set(V46, Ti45, scratch);

    scratch = -screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*state.rho - screened_rates(k_V46_to_He4_Sc42_derived) - screened_rates(k_V46_to_p_Ti45_derived);
    jac.set(V46, V46, scratch);

    jac.set(V46, V47, 0.0);

    jac.set(V46, Cr48, 0.0);

    scratch = screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(H1)*state.rho;
    jac.set(V46, Cr49, scratch);

    scratch = screened_rates(k_Mn50_to_He4_V46_derived);
    jac.set(V46, Mn50, scratch);

    jac.set(V46, Mn51, 0.0);

    jac.set(V46, Mn52, 0.0);

    jac.set(V46, Fe52, 0.0);

    jac.set(V46, Fe53, 0.0);

    jac.set(V46, Fe54, 0.0);

    jac.set(V46, Fe55, 0.0);

    jac.set(V46, Fe56, 0.0);

    jac.set(V46, Co55, 0.0);

    jac.set(V46, Co56, 0.0);

    jac.set(V46, Co57, 0.0);

    jac.set(V46, Ni56, 0.0);

    jac.set(V46, Ni57, 0.0);

    jac.set(V46, Ni58, 0.0);

    jac.set(V46, Cu59, 0.0);

    jac.set(V46, Zn60, 0.0);

    jac.set(V47, N, 0.0);

    scratch = -screened_rates(k_p_V47_to_Cr48_reaclib)*Y(V47)*state.rho - screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(V47)*state.rho;
    jac.set(V47, H1, scratch);

    scratch = screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(Sc43)*state.rho + screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(Ti44)*state.rho - screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(V47)*state.rho;
    jac.set(V47, He4, scratch);

    jac.set(V47, C12, 0.0);

    jac.set(V47, C13, 0.0);

    jac.set(V47, C14, 0.0);

    jac.set(V47, N13, 0.0);

    jac.set(V47, N14, 0.0);

    jac.set(V47, N15, 0.0);

    jac.set(V47, O16, 0.0);

    jac.set(V47, O17, 0.0);

    jac.set(V47, O18, 0.0);

    jac.set(V47, F18, 0.0);

    jac.set(V47, F19, 0.0);

    jac.set(V47, Ne20, 0.0);

    jac.set(V47, Ne21, 0.0);

    jac.set(V47, Ne22, 0.0);

    jac.set(V47, Na22, 0.0);

    jac.set(V47, Na23, 0.0);

    jac.set(V47, Mg24, 0.0);

    jac.set(V47, Mg25, 0.0);

    jac.set(V47, Al26, 0.0);

    jac.set(V47, Al27, 0.0);

    jac.set(V47, Si28, 0.0);

    jac.set(V47, Si29, 0.0);

    jac.set(V47, P30, 0.0);

    jac.set(V47, P31, 0.0);

    jac.set(V47, S32, 0.0);

    jac.set(V47, S33, 0.0);

    jac.set(V47, Cl34, 0.0);

    jac.set(V47, Cl35, 0.0);

    jac.set(V47, Ar36, 0.0);

    jac.set(V47, Ar37, 0.0);

    jac.set(V47, K38, 0.0);

    jac.set(V47, K39, 0.0);

    jac.set(V47, Ca40, 0.0);

    jac.set(V47, Ca41, 0.0);

    jac.set(V47, Sc42, 0.0);

    scratch = screened_rates(k_He4_Sc43_to_V47_reaclib)*Y(He4)*state.rho;
    jac.set(V47, Sc43, scratch);

    scratch = screened_rates(k_He4_Ti44_to_p_V47_reaclib)*Y(He4)*state.rho;
    jac.set(V47, Ti44, scratch);

    jac.set(V47, Ti45, 0.0);

    jac.set(V47, V46, 0.0);

    scratch = -screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(He4)*state.rho - screened_rates(k_V47_to_He4_Sc43_derived) - screened_rates(k_p_V47_to_Cr48_reaclib)*Y(H1)*state.rho - screened_rates(k_p_V47_to_He4_Ti44_derived)*Y(H1)*state.rho;
    jac.set(V47, V47, scratch);

    scratch = screened_rates(k_Cr48_to_p_V47_derived);
    jac.set(V47, Cr48, scratch);

    jac.set(V47, Cr49, 0.0);

    jac.set(V47, Mn50, 0.0);

    scratch = screened_rates(k_Mn51_to_He4_V47_derived);
    jac.set(V47, Mn51, scratch);

    jac.set(V47, Mn52, 0.0);

    jac.set(V47, Fe52, 0.0);

    jac.set(V47, Fe53, 0.0);

    jac.set(V47, Fe54, 0.0);

    jac.set(V47, Fe55, 0.0);

    jac.set(V47, Fe56, 0.0);

    jac.set(V47, Co55, 0.0);

    jac.set(V47, Co56, 0.0);

    jac.set(V47, Co57, 0.0);

    jac.set(V47, Ni56, 0.0);

    jac.set(V47, Ni57, 0.0);

    jac.set(V47, Ni58, 0.0);

    jac.set(V47, Cu59, 0.0);

    jac.set(V47, Zn60, 0.0);

    jac.set(Cr48, N, 0.0);

    scratch = screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho + screened_rates(k_p_V47_to_Cr48_reaclib)*Y(V47)*state.rho;
    jac.set(Cr48, H1, scratch);

    scratch = -screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(Cr48)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*state.rho + screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    jac.set(Cr48, C12, 0.0);

    jac.set(Cr48, C13, 0.0);

    jac.set(Cr48, C14, 0.0);

    jac.set(Cr48, N13, 0.0);

    jac.set(Cr48, N14, 0.0);

    jac.set(Cr48, N15, 0.0);

    jac.set(Cr48, O16, 0.0);

    jac.set(Cr48, O17, 0.0);

    jac.set(Cr48, O18, 0.0);

    jac.set(Cr48, F18, 0.0);

    jac.set(Cr48, F19, 0.0);

    jac.set(Cr48, Ne20, 0.0);

    jac.set(Cr48, Ne21, 0.0);

    jac.set(Cr48, Ne22, 0.0);

    jac.set(Cr48, Na22, 0.0);

    jac.set(Cr48, Na23, 0.0);

    jac.set(Cr48, Mg24, 0.0);

    jac.set(Cr48, Mg25, 0.0);

    jac.set(Cr48, Al26, 0.0);

    jac.set(Cr48, Al27, 0.0);

    jac.set(Cr48, Si28, 0.0);

    jac.set(Cr48, Si29, 0.0);

    jac.set(Cr48, P30, 0.0);

    jac.set(Cr48, P31, 0.0);

    jac.set(Cr48, S32, 0.0);

    jac.set(Cr48, S33, 0.0);

    jac.set(Cr48, Cl34, 0.0);

    jac.set(Cr48, Cl35, 0.0);

    jac.set(Cr48, Ar36, 0.0);

    jac.set(Cr48, Ar37, 0.0);

    jac.set(Cr48, K38, 0.0);

    jac.set(Cr48, K39, 0.0);

    jac.set(Cr48, Ca40, 0.0);

    jac.set(Cr48, Ca41, 0.0);

    jac.set(Cr48, Sc42, 0.0);

    jac.set(Cr48, Sc43, 0.0);

    scratch = screened_rates(k_He4_Ti44_to_Cr48_reaclib)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    jac.set(Cr48, Ti45, 0.0);

    jac.set(Cr48, V46, 0.0);

    scratch = screened_rates(k_p_V47_to_Cr48_reaclib)*Y(H1)*state.rho;
    jac.set(Cr48, V47, scratch);

    scratch = -screened_rates(k_Cr48_to_He4_Ti44_derived) - screened_rates(k_Cr48_to_p_V47_derived) - screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(He4)*state.rho;
    jac.set(Cr48, Cr48, scratch);

    jac.set(Cr48, Cr49, 0.0);

    jac.set(Cr48, Mn50, 0.0);

    scratch = screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(H1)*state.rho;
    jac.set(Cr48, Mn51, scratch);

    jac.set(Cr48, Mn52, 0.0);

    scratch = screened_rates(k_Fe52_to_He4_Cr48_derived);
    jac.set(Cr48, Fe52, scratch);

    jac.set(Cr48, Fe53, 0.0);

    jac.set(Cr48, Fe54, 0.0);

    jac.set(Cr48, Fe55, 0.0);

    jac.set(Cr48, Fe56, 0.0);

    jac.set(Cr48, Co55, 0.0);

    jac.set(Cr48, Co56, 0.0);

    jac.set(Cr48, Co57, 0.0);

    jac.set(Cr48, Ni56, 0.0);

    jac.set(Cr48, Ni57, 0.0);

    jac.set(Cr48, Ni58, 0.0);

    jac.set(Cr48, Cu59, 0.0);

    jac.set(Cr48, Zn60, 0.0);

    jac.set(Cr49, N, 0.0);

    scratch = -screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(Cr49)*state.rho - screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(Cr49)*state.rho + screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*state.rho;
    jac.set(Cr49, H1, scratch);

    scratch = -screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*state.rho + screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(Ti45)*state.rho + screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(V46)*state.rho;
    jac.set(Cr49, He4, scratch);

    jac.set(Cr49, C12, 0.0);

    jac.set(Cr49, C13, 0.0);

    jac.set(Cr49, C14, 0.0);

    jac.set(Cr49, N13, 0.0);

    jac.set(Cr49, N14, 0.0);

    jac.set(Cr49, N15, 0.0);

    jac.set(Cr49, O16, 0.0);

    jac.set(Cr49, O17, 0.0);

    jac.set(Cr49, O18, 0.0);

    jac.set(Cr49, F18, 0.0);

    jac.set(Cr49, F19, 0.0);

    jac.set(Cr49, Ne20, 0.0);

    jac.set(Cr49, Ne21, 0.0);

    jac.set(Cr49, Ne22, 0.0);

    jac.set(Cr49, Na22, 0.0);

    jac.set(Cr49, Na23, 0.0);

    jac.set(Cr49, Mg24, 0.0);

    jac.set(Cr49, Mg25, 0.0);

    jac.set(Cr49, Al26, 0.0);

    jac.set(Cr49, Al27, 0.0);

    jac.set(Cr49, Si28, 0.0);

    jac.set(Cr49, Si29, 0.0);

    jac.set(Cr49, P30, 0.0);

    jac.set(Cr49, P31, 0.0);

    jac.set(Cr49, S32, 0.0);

    jac.set(Cr49, S33, 0.0);

    jac.set(Cr49, Cl34, 0.0);

    jac.set(Cr49, Cl35, 0.0);

    jac.set(Cr49, Ar36, 0.0);

    jac.set(Cr49, Ar37, 0.0);

    jac.set(Cr49, K38, 0.0);

    jac.set(Cr49, K39, 0.0);

    jac.set(Cr49, Ca40, 0.0);

    jac.set(Cr49, Ca41, 0.0);

    jac.set(Cr49, Sc42, 0.0);

    jac.set(Cr49, Sc43, 0.0);

    jac.set(Cr49, Ti44, 0.0);

    scratch = screened_rates(k_He4_Ti45_to_Cr49_reaclib)*Y(He4)*state.rho;
    jac.set(Cr49, Ti45, scratch);

    scratch = screened_rates(k_He4_V46_to_p_Cr49_reaclib)*Y(He4)*state.rho;
    jac.set(Cr49, V46, scratch);

    jac.set(Cr49, V47, 0.0);

    jac.set(Cr49, Cr48, 0.0);

    scratch = -screened_rates(k_Cr49_to_He4_Ti45_derived) - screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(He4)*state.rho - screened_rates(k_p_Cr49_to_He4_V46_derived)*Y(H1)*state.rho - screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(H1)*state.rho;
    jac.set(Cr49, Cr49, scratch);

    scratch = screened_rates(k_Mn50_to_p_Cr49_derived);
    jac.set(Cr49, Mn50, scratch);

    jac.set(Cr49, Mn51, 0.0);

    scratch = screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(H1)*state.rho;
    jac.set(Cr49, Mn52, scratch);

    jac.set(Cr49, Fe52, 0.0);

    jac.set(Cr49, Fe53, 0.0);

    jac.set(Cr49, Fe54, 0.0);

    jac.set(Cr49, Fe55, 0.0);

    jac.set(Cr49, Fe56, 0.0);

    jac.set(Cr49, Co55, 0.0);

    jac.set(Cr49, Co56, 0.0);

    jac.set(Cr49, Co57, 0.0);

    jac.set(Cr49, Ni56, 0.0);

    jac.set(Cr49, Ni57, 0.0);

    jac.set(Cr49, Ni58, 0.0);

    jac.set(Cr49, Cu59, 0.0);

    jac.set(Cr49, Zn60, 0.0);

    scratch = -screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(Mn50)*state.rho;
    jac.set(Mn50, N, scratch);

    scratch = screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(Cr49)*state.rho + screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*state.rho;
    jac.set(Mn50, H1, scratch);

    scratch = -screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(Mn50)*state.rho + screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(V46)*state.rho;
    jac.set(Mn50, He4, scratch);

    jac.set(Mn50, C12, 0.0);

    jac.set(Mn50, C13, 0.0);

    jac.set(Mn50, C14, 0.0);

    jac.set(Mn50, N13, 0.0);

    jac.set(Mn50, N14, 0.0);

    jac.set(Mn50, N15, 0.0);

    jac.set(Mn50, O16, 0.0);

    jac.set(Mn50, O17, 0.0);

    jac.set(Mn50, O18, 0.0);

    jac.set(Mn50, F18, 0.0);

    jac.set(Mn50, F19, 0.0);

    jac.set(Mn50, Ne20, 0.0);

    jac.set(Mn50, Ne21, 0.0);

    jac.set(Mn50, Ne22, 0.0);

    jac.set(Mn50, Na22, 0.0);

    jac.set(Mn50, Na23, 0.0);

    jac.set(Mn50, Mg24, 0.0);

    jac.set(Mn50, Mg25, 0.0);

    jac.set(Mn50, Al26, 0.0);

    jac.set(Mn50, Al27, 0.0);

    jac.set(Mn50, Si28, 0.0);

    jac.set(Mn50, Si29, 0.0);

    jac.set(Mn50, P30, 0.0);

    jac.set(Mn50, P31, 0.0);

    jac.set(Mn50, S32, 0.0);

    jac.set(Mn50, S33, 0.0);

    jac.set(Mn50, Cl34, 0.0);

    jac.set(Mn50, Cl35, 0.0);

    jac.set(Mn50, Ar36, 0.0);

    jac.set(Mn50, Ar37, 0.0);

    jac.set(Mn50, K38, 0.0);

    jac.set(Mn50, K39, 0.0);

    jac.set(Mn50, Ca40, 0.0);

    jac.set(Mn50, Ca41, 0.0);

    jac.set(Mn50, Sc42, 0.0);

    jac.set(Mn50, Sc43, 0.0);

    jac.set(Mn50, Ti44, 0.0);

    jac.set(Mn50, Ti45, 0.0);

    scratch = screened_rates(k_He4_V46_to_Mn50_reaclib)*Y(He4)*state.rho;
    jac.set(Mn50, V46, scratch);

    jac.set(Mn50, V47, 0.0);

    jac.set(Mn50, Cr48, 0.0);

    scratch = screened_rates(k_p_Cr49_to_Mn50_reaclib)*Y(H1)*state.rho;
    jac.set(Mn50, Cr49, scratch);

    scratch = -screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*state.rho - screened_rates(k_Mn50_to_He4_V46_derived) - screened_rates(k_Mn50_to_p_Cr49_derived) - screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(N)*state.rho;
    jac.set(Mn50, Mn50, scratch);

    scratch = screened_rates(k_Mn51_to_n_Mn50_derived);
    jac.set(Mn50, Mn51, scratch);

    jac.set(Mn50, Mn52, 0.0);

    jac.set(Mn50, Fe52, 0.0);

    scratch = screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(H1)*state.rho;
    jac.set(Mn50, Fe53, scratch);

    jac.set(Mn50, Fe54, 0.0);

    jac.set(Mn50, Fe55, 0.0);

    jac.set(Mn50, Fe56, 0.0);

    jac.set(Mn50, Co55, 0.0);

    jac.set(Mn50, Co56, 0.0);

    jac.set(Mn50, Co57, 0.0);

    jac.set(Mn50, Ni56, 0.0);

    jac.set(Mn50, Ni57, 0.0);

    jac.set(Mn50, Ni58, 0.0);

    jac.set(Mn50, Cu59, 0.0);

    jac.set(Mn50, Zn60, 0.0);

    scratch = screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(Mn50)*state.rho - screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(Mn51)*state.rho;
    jac.set(Mn51, N, scratch);

    scratch = screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho - screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(Mn51)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho;
    jac.set(Mn51, H1, scratch);

    scratch = screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(Cr48)*state.rho - screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(Mn51)*state.rho - screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(Mn51)*state.rho + screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(V47)*state.rho;
    jac.set(Mn51, He4, scratch);

    jac.set(Mn51, C12, 0.0);

    jac.set(Mn51, C13, 0.0);

    jac.set(Mn51, C14, 0.0);

    jac.set(Mn51, N13, 0.0);

    jac.set(Mn51, N14, 0.0);

    jac.set(Mn51, N15, 0.0);

    jac.set(Mn51, O16, 0.0);

    jac.set(Mn51, O17, 0.0);

    jac.set(Mn51, O18, 0.0);

    jac.set(Mn51, F18, 0.0);

    jac.set(Mn51, F19, 0.0);

    jac.set(Mn51, Ne20, 0.0);

    jac.set(Mn51, Ne21, 0.0);

    jac.set(Mn51, Ne22, 0.0);

    jac.set(Mn51, Na22, 0.0);

    jac.set(Mn51, Na23, 0.0);

    jac.set(Mn51, Mg24, 0.0);

    jac.set(Mn51, Mg25, 0.0);

    jac.set(Mn51, Al26, 0.0);

    jac.set(Mn51, Al27, 0.0);

    jac.set(Mn51, Si28, 0.0);

    jac.set(Mn51, Si29, 0.0);

    jac.set(Mn51, P30, 0.0);

    jac.set(Mn51, P31, 0.0);

    jac.set(Mn51, S32, 0.0);

    jac.set(Mn51, S33, 0.0);

    jac.set(Mn51, Cl34, 0.0);

    jac.set(Mn51, Cl35, 0.0);

    jac.set(Mn51, Ar36, 0.0);

    jac.set(Mn51, Ar37, 0.0);

    jac.set(Mn51, K38, 0.0);

    jac.set(Mn51, K39, 0.0);

    jac.set(Mn51, Ca40, 0.0);

    jac.set(Mn51, Ca41, 0.0);

    jac.set(Mn51, Sc42, 0.0);

    jac.set(Mn51, Sc43, 0.0);

    jac.set(Mn51, Ti44, 0.0);

    jac.set(Mn51, Ti45, 0.0);

    jac.set(Mn51, V46, 0.0);

    scratch = screened_rates(k_He4_V47_to_Mn51_reaclib)*Y(He4)*state.rho;
    jac.set(Mn51, V47, scratch);

    scratch = screened_rates(k_He4_Cr48_to_p_Mn51_reaclib)*Y(He4)*state.rho;
    jac.set(Mn51, Cr48, scratch);

    jac.set(Mn51, Cr49, 0.0);

    scratch = screened_rates(k_n_Mn50_to_Mn51_reaclib)*Y(N)*state.rho;
    jac.set(Mn51, Mn50, scratch);

    scratch = -screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*state.rho - screened_rates(k_Mn51_to_He4_V47_derived) - screened_rates(k_Mn51_to_n_Mn50_derived) - screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(N)*state.rho - screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48_derived)*Y(H1)*state.rho;
    jac.set(Mn51, Mn51, scratch);

    scratch = screened_rates(k_Mn52_to_n_Mn51_derived);
    jac.set(Mn51, Mn52, scratch);

    scratch = screened_rates(k_Fe52_to_p_Mn51_derived);
    jac.set(Mn51, Fe52, scratch);

    jac.set(Mn51, Fe53, 0.0);

    scratch = screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(H1)*state.rho;
    jac.set(Mn51, Fe54, scratch);

    jac.set(Mn51, Fe55, 0.0);

    jac.set(Mn51, Fe56, 0.0);

    scratch = screened_rates(k_Co55_to_He4_Mn51_derived);
    jac.set(Mn51, Co55, scratch);

    jac.set(Mn51, Co56, 0.0);

    jac.set(Mn51, Co57, 0.0);

    jac.set(Mn51, Ni56, 0.0);

    jac.set(Mn51, Ni57, 0.0);

    jac.set(Mn51, Ni58, 0.0);

    jac.set(Mn51, Cu59, 0.0);

    jac.set(Mn51, Zn60, 0.0);

    scratch = screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*state.rho + screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*state.rho + screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(Mn51)*state.rho;
    jac.set(Mn52, N, scratch);

    scratch = screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*state.rho - screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(Mn52)*state.rho - screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(Mn52)*state.rho - screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*state.rho;
    jac.set(Mn52, H1, scratch);

    scratch = screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(Cr49)*state.rho - screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(Mn52)*state.rho - screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(Mn52)*state.rho - screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(Mn52)*state.rho;
    jac.set(Mn52, He4, scratch);

    jac.set(Mn52, C12, 0.0);

    jac.set(Mn52, C13, 0.0);

    jac.set(Mn52, C14, 0.0);

    jac.set(Mn52, N13, 0.0);

    jac.set(Mn52, N14, 0.0);

    jac.set(Mn52, N15, 0.0);

    jac.set(Mn52, O16, 0.0);

    jac.set(Mn52, O17, 0.0);

    jac.set(Mn52, O18, 0.0);

    jac.set(Mn52, F18, 0.0);

    jac.set(Mn52, F19, 0.0);

    jac.set(Mn52, Ne20, 0.0);

    jac.set(Mn52, Ne21, 0.0);

    jac.set(Mn52, Ne22, 0.0);

    jac.set(Mn52, Na22, 0.0);

    jac.set(Mn52, Na23, 0.0);

    jac.set(Mn52, Mg24, 0.0);

    jac.set(Mn52, Mg25, 0.0);

    jac.set(Mn52, Al26, 0.0);

    jac.set(Mn52, Al27, 0.0);

    jac.set(Mn52, Si28, 0.0);

    jac.set(Mn52, Si29, 0.0);

    jac.set(Mn52, P30, 0.0);

    jac.set(Mn52, P31, 0.0);

    jac.set(Mn52, S32, 0.0);

    jac.set(Mn52, S33, 0.0);

    jac.set(Mn52, Cl34, 0.0);

    jac.set(Mn52, Cl35, 0.0);

    jac.set(Mn52, Ar36, 0.0);

    jac.set(Mn52, Ar37, 0.0);

    jac.set(Mn52, K38, 0.0);

    jac.set(Mn52, K39, 0.0);

    jac.set(Mn52, Ca40, 0.0);

    jac.set(Mn52, Ca41, 0.0);

    jac.set(Mn52, Sc42, 0.0);

    jac.set(Mn52, Sc43, 0.0);

    jac.set(Mn52, Ti44, 0.0);

    jac.set(Mn52, Ti45, 0.0);

    jac.set(Mn52, V46, 0.0);

    jac.set(Mn52, V47, 0.0);

    jac.set(Mn52, Cr48, 0.0);

    scratch = screened_rates(k_He4_Cr49_to_p_Mn52_reaclib)*Y(He4)*state.rho;
    jac.set(Mn52, Cr49, scratch);

    jac.set(Mn52, Mn50, 0.0);

    scratch = screened_rates(k_n_Mn51_to_Mn52_reaclib)*Y(N)*state.rho;
    jac.set(Mn52, Mn51, scratch);

    scratch = -screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*state.rho - screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*state.rho - screened_rates(k_Mn52_to_Fe52_weaktab) - screened_rates(k_Mn52_to_n_Mn51_derived) - screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Mn52_to_He4_Cr49_derived)*Y(H1)*state.rho - screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(H1)*state.rho;
    jac.set(Mn52, Mn52, scratch);

    scratch = screened_rates(k_Fe52_to_Mn52_weaktab) + screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(N)*state.rho;
    jac.set(Mn52, Fe52, scratch);

    scratch = screened_rates(k_Fe53_to_p_Mn52_derived);
    jac.set(Mn52, Fe53, scratch);

    jac.set(Mn52, Fe54, 0.0);

    scratch = screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(H1)*state.rho;
    jac.set(Mn52, Fe55, scratch);

    jac.set(Mn52, Fe56, 0.0);

    scratch = screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(N)*state.rho;
    jac.set(Mn52, Co55, scratch);

    scratch = screened_rates(k_Co56_to_He4_Mn52_derived);
    jac.set(Mn52, Co56, scratch);

    jac.set(Mn52, Co57, 0.0);

    jac.set(Mn52, Ni56, 0.0);

    jac.set(Mn52, Ni57, 0.0);

    jac.set(Mn52, Ni58, 0.0);

    jac.set(Mn52, Cu59, 0.0);

    jac.set(Mn52, Zn60, 0.0);

    scratch = -screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(Fe52)*state.rho - screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(Fe52)*state.rho;
    jac.set(Fe52, N, scratch);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho + screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(Mn51)*state.rho + screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(Mn52)*state.rho;
    jac.set(Fe52, H1, scratch);

    scratch = screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(Cr48)*state.rho - screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    jac.set(Fe52, C12, 0.0);

    jac.set(Fe52, C13, 0.0);

    jac.set(Fe52, C14, 0.0);

    jac.set(Fe52, N13, 0.0);

    jac.set(Fe52, N14, 0.0);

    jac.set(Fe52, N15, 0.0);

    jac.set(Fe52, O16, 0.0);

    jac.set(Fe52, O17, 0.0);

    jac.set(Fe52, O18, 0.0);

    jac.set(Fe52, F18, 0.0);

    jac.set(Fe52, F19, 0.0);

    jac.set(Fe52, Ne20, 0.0);

    jac.set(Fe52, Ne21, 0.0);

    jac.set(Fe52, Ne22, 0.0);

    jac.set(Fe52, Na22, 0.0);

    jac.set(Fe52, Na23, 0.0);

    jac.set(Fe52, Mg24, 0.0);

    jac.set(Fe52, Mg25, 0.0);

    jac.set(Fe52, Al26, 0.0);

    jac.set(Fe52, Al27, 0.0);

    jac.set(Fe52, Si28, 0.0);

    jac.set(Fe52, Si29, 0.0);

    jac.set(Fe52, P30, 0.0);

    jac.set(Fe52, P31, 0.0);

    jac.set(Fe52, S32, 0.0);

    jac.set(Fe52, S33, 0.0);

    jac.set(Fe52, Cl34, 0.0);

    jac.set(Fe52, Cl35, 0.0);

    jac.set(Fe52, Ar36, 0.0);

    jac.set(Fe52, Ar37, 0.0);

    jac.set(Fe52, K38, 0.0);

    jac.set(Fe52, K39, 0.0);

    jac.set(Fe52, Ca40, 0.0);

    jac.set(Fe52, Ca41, 0.0);

    jac.set(Fe52, Sc42, 0.0);

    jac.set(Fe52, Sc43, 0.0);

    jac.set(Fe52, Ti44, 0.0);

    jac.set(Fe52, Ti45, 0.0);

    jac.set(Fe52, V46, 0.0);

    jac.set(Fe52, V47, 0.0);

    scratch = screened_rates(k_He4_Cr48_to_Fe52_reaclib)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    jac.set(Fe52, Cr49, 0.0);

    jac.set(Fe52, Mn50, 0.0);

    scratch = screened_rates(k_p_Mn51_to_Fe52_reaclib)*Y(H1)*state.rho;
    jac.set(Fe52, Mn51, scratch);

    scratch = screened_rates(k_Mn52_to_Fe52_weaktab) + screened_rates(k_p_Mn52_to_n_Fe52_derived)*Y(H1)*state.rho;
    jac.set(Fe52, Mn52, scratch);

    scratch = -screened_rates(k_Fe52_to_He4_Cr48_derived) - screened_rates(k_Fe52_to_Mn52_weaktab) - screened_rates(k_Fe52_to_p_Mn51_derived) - screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(He4)*state.rho - screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(N)*state.rho - screened_rates(k_n_Fe52_to_p_Mn52_reaclib)*Y(N)*state.rho;
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_Fe53_to_n_Fe52_derived);
    jac.set(Fe52, Fe53, scratch);

    jac.set(Fe52, Fe54, 0.0);

    jac.set(Fe52, Fe55, 0.0);

    jac.set(Fe52, Fe56, 0.0);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho;
    jac.set(Fe52, Co55, scratch);

    jac.set(Fe52, Co56, 0.0);

    jac.set(Fe52, Co57, 0.0);

    scratch = screened_rates(k_Ni56_to_He4_Fe52_derived);
    jac.set(Fe52, Ni56, scratch);

    jac.set(Fe52, Ni57, 0.0);

    jac.set(Fe52, Ni58, 0.0);

    jac.set(Fe52, Cu59, 0.0);

    jac.set(Fe52, Zn60, 0.0);

    scratch = screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(Fe52)*state.rho - screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(Fe53)*state.rho + screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*state.rho;
    jac.set(Fe53, N, scratch);

    scratch = screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*state.rho - screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(Fe53)*state.rho + screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(Mn52)*state.rho;
    jac.set(Fe53, H1, scratch);

    scratch = -screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(Fe53)*state.rho - screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*state.rho - screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*state.rho + screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(Mn50)*state.rho;
    jac.set(Fe53, He4, scratch);

    jac.set(Fe53, C12, 0.0);

    jac.set(Fe53, C13, 0.0);

    jac.set(Fe53, C14, 0.0);

    jac.set(Fe53, N13, 0.0);

    jac.set(Fe53, N14, 0.0);

    jac.set(Fe53, N15, 0.0);

    jac.set(Fe53, O16, 0.0);

    jac.set(Fe53, O17, 0.0);

    jac.set(Fe53, O18, 0.0);

    jac.set(Fe53, F18, 0.0);

    jac.set(Fe53, F19, 0.0);

    jac.set(Fe53, Ne20, 0.0);

    jac.set(Fe53, Ne21, 0.0);

    jac.set(Fe53, Ne22, 0.0);

    jac.set(Fe53, Na22, 0.0);

    jac.set(Fe53, Na23, 0.0);

    jac.set(Fe53, Mg24, 0.0);

    jac.set(Fe53, Mg25, 0.0);

    jac.set(Fe53, Al26, 0.0);

    jac.set(Fe53, Al27, 0.0);

    jac.set(Fe53, Si28, 0.0);

    jac.set(Fe53, Si29, 0.0);

    jac.set(Fe53, P30, 0.0);

    jac.set(Fe53, P31, 0.0);

    jac.set(Fe53, S32, 0.0);

    jac.set(Fe53, S33, 0.0);

    jac.set(Fe53, Cl34, 0.0);

    jac.set(Fe53, Cl35, 0.0);

    jac.set(Fe53, Ar36, 0.0);

    jac.set(Fe53, Ar37, 0.0);

    jac.set(Fe53, K38, 0.0);

    jac.set(Fe53, K39, 0.0);

    jac.set(Fe53, Ca40, 0.0);

    jac.set(Fe53, Ca41, 0.0);

    jac.set(Fe53, Sc42, 0.0);

    jac.set(Fe53, Sc43, 0.0);

    jac.set(Fe53, Ti44, 0.0);

    jac.set(Fe53, Ti45, 0.0);

    jac.set(Fe53, V46, 0.0);

    jac.set(Fe53, V47, 0.0);

    jac.set(Fe53, Cr48, 0.0);

    jac.set(Fe53, Cr49, 0.0);

    scratch = screened_rates(k_He4_Mn50_to_p_Fe53_reaclib)*Y(He4)*state.rho;
    jac.set(Fe53, Mn50, scratch);

    jac.set(Fe53, Mn51, 0.0);

    scratch = screened_rates(k_p_Mn52_to_Fe53_reaclib)*Y(H1)*state.rho;
    jac.set(Fe53, Mn52, scratch);

    scratch = screened_rates(k_n_Fe52_to_Fe53_reaclib)*Y(N)*state.rho;
    jac.set(Fe53, Fe52, scratch);

    scratch = -screened_rates(k_Fe53_to_n_Fe52_derived) - screened_rates(k_Fe53_to_p_Mn52_derived) - screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(He4)*state.rho - screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(He4)*state.rho - screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(N)*state.rho - screened_rates(k_p_Fe53_to_He4_Mn50_derived)*Y(H1)*state.rho;
    jac.set(Fe53, Fe53, scratch);

    scratch = screened_rates(k_Fe54_to_n_Fe53_derived);
    jac.set(Fe53, Fe54, scratch);

    jac.set(Fe53, Fe55, 0.0);

    jac.set(Fe53, Fe56, 0.0);

    jac.set(Fe53, Co55, 0.0);

    scratch = screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(H1)*state.rho;
    jac.set(Fe53, Co56, scratch);

    jac.set(Fe53, Co57, 0.0);

    scratch = screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(N)*state.rho;
    jac.set(Fe53, Ni56, scratch);

    scratch = screened_rates(k_Ni57_to_He4_Fe53_derived);
    jac.set(Fe53, Ni57, scratch);

    jac.set(Fe53, Ni58, 0.0);

    jac.set(Fe53, Cu59, 0.0);

    jac.set(Fe53, Zn60, 0.0);

    scratch = screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(Fe53)*state.rho - screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(Fe54)*state.rho + screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*state.rho;
    jac.set(Fe54, N, scratch);

    scratch = screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*state.rho - screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(Fe54)*state.rho - screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho;
    jac.set(Fe54, H1, scratch);

    scratch = -screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(Fe54)*state.rho - screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*state.rho - screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*state.rho + screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(Mn51)*state.rho;
    jac.set(Fe54, He4, scratch);

    jac.set(Fe54, C12, 0.0);

    jac.set(Fe54, C13, 0.0);

    jac.set(Fe54, C14, 0.0);

    jac.set(Fe54, N13, 0.0);

    jac.set(Fe54, N14, 0.0);

    jac.set(Fe54, N15, 0.0);

    jac.set(Fe54, O16, 0.0);

    jac.set(Fe54, O17, 0.0);

    jac.set(Fe54, O18, 0.0);

    jac.set(Fe54, F18, 0.0);

    jac.set(Fe54, F19, 0.0);

    jac.set(Fe54, Ne20, 0.0);

    jac.set(Fe54, Ne21, 0.0);

    jac.set(Fe54, Ne22, 0.0);

    jac.set(Fe54, Na22, 0.0);

    jac.set(Fe54, Na23, 0.0);

    jac.set(Fe54, Mg24, 0.0);

    jac.set(Fe54, Mg25, 0.0);

    jac.set(Fe54, Al26, 0.0);

    jac.set(Fe54, Al27, 0.0);

    jac.set(Fe54, Si28, 0.0);

    jac.set(Fe54, Si29, 0.0);

    jac.set(Fe54, P30, 0.0);

    jac.set(Fe54, P31, 0.0);

    jac.set(Fe54, S32, 0.0);

    jac.set(Fe54, S33, 0.0);

    jac.set(Fe54, Cl34, 0.0);

    jac.set(Fe54, Cl35, 0.0);

    jac.set(Fe54, Ar36, 0.0);

    jac.set(Fe54, Ar37, 0.0);

    jac.set(Fe54, K38, 0.0);

    jac.set(Fe54, K39, 0.0);

    jac.set(Fe54, Ca40, 0.0);

    jac.set(Fe54, Ca41, 0.0);

    jac.set(Fe54, Sc42, 0.0);

    jac.set(Fe54, Sc43, 0.0);

    jac.set(Fe54, Ti44, 0.0);

    jac.set(Fe54, Ti45, 0.0);

    jac.set(Fe54, V46, 0.0);

    jac.set(Fe54, V47, 0.0);

    jac.set(Fe54, Cr48, 0.0);

    jac.set(Fe54, Cr49, 0.0);

    jac.set(Fe54, Mn50, 0.0);

    scratch = screened_rates(k_He4_Mn51_to_p_Fe54_reaclib)*Y(He4)*state.rho;
    jac.set(Fe54, Mn51, scratch);

    jac.set(Fe54, Mn52, 0.0);

    jac.set(Fe54, Fe52, 0.0);

    scratch = screened_rates(k_n_Fe53_to_Fe54_reaclib)*Y(N)*state.rho;
    jac.set(Fe54, Fe53, scratch);

    scratch = -screened_rates(k_Fe54_to_n_Fe53_derived) - screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(He4)*state.rho - screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(He4)*state.rho - screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(N)*state.rho - screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Fe54_to_He4_Mn51_derived)*Y(H1)*state.rho;
    jac.set(Fe54, Fe54, scratch);

    scratch = screened_rates(k_Fe55_to_n_Fe54_derived);
    jac.set(Fe54, Fe55, scratch);

    jac.set(Fe54, Fe56, 0.0);

    scratch = screened_rates(k_Co55_to_p_Fe54_derived);
    jac.set(Fe54, Co55, scratch);

    jac.set(Fe54, Co56, 0.0);

    scratch = screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(H1)*state.rho;
    jac.set(Fe54, Co57, scratch);

    jac.set(Fe54, Ni56, 0.0);

    scratch = screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(N)*state.rho;
    jac.set(Fe54, Ni57, scratch);

    scratch = screened_rates(k_Ni58_to_He4_Fe54_derived);
    jac.set(Fe54, Ni58, scratch);

    jac.set(Fe54, Cu59, 0.0);

    jac.set(Fe54, Zn60, 0.0);

    scratch = screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*state.rho + screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(Fe54)*state.rho - screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(Fe55)*state.rho + screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*state.rho;
    jac.set(Fe55, N, scratch);

    scratch = -screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(Fe55)*state.rho - screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(Fe55)*state.rho - screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*state.rho;
    jac.set(Fe55, H1, scratch);

    scratch = -screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*state.rho + screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(Mn52)*state.rho;
    jac.set(Fe55, He4, scratch);

    jac.set(Fe55, C12, 0.0);

    jac.set(Fe55, C13, 0.0);

    jac.set(Fe55, C14, 0.0);

    jac.set(Fe55, N13, 0.0);

    jac.set(Fe55, N14, 0.0);

    jac.set(Fe55, N15, 0.0);

    jac.set(Fe55, O16, 0.0);

    jac.set(Fe55, O17, 0.0);

    jac.set(Fe55, O18, 0.0);

    jac.set(Fe55, F18, 0.0);

    jac.set(Fe55, F19, 0.0);

    jac.set(Fe55, Ne20, 0.0);

    jac.set(Fe55, Ne21, 0.0);

    jac.set(Fe55, Ne22, 0.0);

    jac.set(Fe55, Na22, 0.0);

    jac.set(Fe55, Na23, 0.0);

    jac.set(Fe55, Mg24, 0.0);

    jac.set(Fe55, Mg25, 0.0);

    jac.set(Fe55, Al26, 0.0);

    jac.set(Fe55, Al27, 0.0);

    jac.set(Fe55, Si28, 0.0);

    jac.set(Fe55, Si29, 0.0);

    jac.set(Fe55, P30, 0.0);

    jac.set(Fe55, P31, 0.0);

    jac.set(Fe55, S32, 0.0);

    jac.set(Fe55, S33, 0.0);

    jac.set(Fe55, Cl34, 0.0);

    jac.set(Fe55, Cl35, 0.0);

    jac.set(Fe55, Ar36, 0.0);

    jac.set(Fe55, Ar37, 0.0);

    jac.set(Fe55, K38, 0.0);

    jac.set(Fe55, K39, 0.0);

    jac.set(Fe55, Ca40, 0.0);

    jac.set(Fe55, Ca41, 0.0);

    jac.set(Fe55, Sc42, 0.0);

    jac.set(Fe55, Sc43, 0.0);

    jac.set(Fe55, Ti44, 0.0);

    jac.set(Fe55, Ti45, 0.0);

    jac.set(Fe55, V46, 0.0);

    jac.set(Fe55, V47, 0.0);

    jac.set(Fe55, Cr48, 0.0);

    jac.set(Fe55, Cr49, 0.0);

    jac.set(Fe55, Mn50, 0.0);

    jac.set(Fe55, Mn51, 0.0);

    scratch = screened_rates(k_He4_Mn52_to_p_Fe55_reaclib)*Y(He4)*state.rho;
    jac.set(Fe55, Mn52, scratch);

    jac.set(Fe55, Fe52, 0.0);

    jac.set(Fe55, Fe53, 0.0);

    scratch = screened_rates(k_n_Fe54_to_Fe55_reaclib)*Y(N)*state.rho;
    jac.set(Fe55, Fe54, scratch);

    scratch = -screened_rates(k_Fe55_to_Co55_weaktab) - screened_rates(k_Fe55_to_n_Fe54_derived) - screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(He4)*state.rho - screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(N)*state.rho - screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Fe55_to_He4_Mn52_derived)*Y(H1)*state.rho - screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(H1)*state.rho;
    jac.set(Fe55, Fe55, scratch);

    scratch = screened_rates(k_Fe56_to_n_Fe55_derived);
    jac.set(Fe55, Fe56, scratch);

    scratch = screened_rates(k_Co55_to_Fe55_weaktab) + screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(N)*state.rho;
    jac.set(Fe55, Co55, scratch);

    scratch = screened_rates(k_Co56_to_p_Fe55_derived);
    jac.set(Fe55, Co56, scratch);

    jac.set(Fe55, Co57, 0.0);

    jac.set(Fe55, Ni56, 0.0);

    jac.set(Fe55, Ni57, 0.0);

    scratch = screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(N)*state.rho;
    jac.set(Fe55, Ni58, scratch);

    jac.set(Fe55, Cu59, 0.0);

    jac.set(Fe55, Zn60, 0.0);

    scratch = screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*state.rho + screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(Fe55)*state.rho;
    jac.set(Fe56, N, scratch);

    scratch = -screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(Fe56)*state.rho - screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho;
    jac.set(Fe56, H1, scratch);

    jac.set(Fe56, He4, 0.0);

    jac.set(Fe56, C12, 0.0);

    jac.set(Fe56, C13, 0.0);

    jac.set(Fe56, C14, 0.0);

    jac.set(Fe56, N13, 0.0);

    jac.set(Fe56, N14, 0.0);

    jac.set(Fe56, N15, 0.0);

    jac.set(Fe56, O16, 0.0);

    jac.set(Fe56, O17, 0.0);

    jac.set(Fe56, O18, 0.0);

    jac.set(Fe56, F18, 0.0);

    jac.set(Fe56, F19, 0.0);

    jac.set(Fe56, Ne20, 0.0);

    jac.set(Fe56, Ne21, 0.0);

    jac.set(Fe56, Ne22, 0.0);

    jac.set(Fe56, Na22, 0.0);

    jac.set(Fe56, Na23, 0.0);

    jac.set(Fe56, Mg24, 0.0);

    jac.set(Fe56, Mg25, 0.0);

    jac.set(Fe56, Al26, 0.0);

    jac.set(Fe56, Al27, 0.0);

    jac.set(Fe56, Si28, 0.0);

    jac.set(Fe56, Si29, 0.0);

    jac.set(Fe56, P30, 0.0);

    jac.set(Fe56, P31, 0.0);

    jac.set(Fe56, S32, 0.0);

    jac.set(Fe56, S33, 0.0);

    jac.set(Fe56, Cl34, 0.0);

    jac.set(Fe56, Cl35, 0.0);

    jac.set(Fe56, Ar36, 0.0);

    jac.set(Fe56, Ar37, 0.0);

    jac.set(Fe56, K38, 0.0);

    jac.set(Fe56, K39, 0.0);

    jac.set(Fe56, Ca40, 0.0);

    jac.set(Fe56, Ca41, 0.0);

    jac.set(Fe56, Sc42, 0.0);

    jac.set(Fe56, Sc43, 0.0);

    jac.set(Fe56, Ti44, 0.0);

    jac.set(Fe56, Ti45, 0.0);

    jac.set(Fe56, V46, 0.0);

    jac.set(Fe56, V47, 0.0);

    jac.set(Fe56, Cr48, 0.0);

    jac.set(Fe56, Cr49, 0.0);

    jac.set(Fe56, Mn50, 0.0);

    jac.set(Fe56, Mn51, 0.0);

    jac.set(Fe56, Mn52, 0.0);

    jac.set(Fe56, Fe52, 0.0);

    jac.set(Fe56, Fe53, 0.0);

    jac.set(Fe56, Fe54, 0.0);

    scratch = screened_rates(k_n_Fe55_to_Fe56_reaclib)*Y(N)*state.rho;
    jac.set(Fe56, Fe55, scratch);

    scratch = -screened_rates(k_Fe56_to_Co56_weaktab) - screened_rates(k_Fe56_to_n_Fe55_derived) - screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(H1)*state.rho;
    jac.set(Fe56, Fe56, scratch);

    jac.set(Fe56, Co55, 0.0);

    scratch = screened_rates(k_Co56_to_Fe56_weaktab) + screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(N)*state.rho;
    jac.set(Fe56, Co56, scratch);

    scratch = screened_rates(k_Co57_to_p_Fe56_derived);
    jac.set(Fe56, Co57, scratch);

    jac.set(Fe56, Ni56, 0.0);

    jac.set(Fe56, Ni57, 0.0);

    jac.set(Fe56, Ni58, 0.0);

    jac.set(Fe56, Cu59, 0.0);

    jac.set(Fe56, Zn60, 0.0);

    scratch = -screened_rates(k_n_Co55_to_Co56_reaclib)*Y(Co55)*state.rho - screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(Co55)*state.rho - screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(Co55)*state.rho;
    jac.set(Co55, N, scratch);

    scratch = -screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(Co55)*state.rho + screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(Fe54)*state.rho + screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(Fe55)*state.rho + screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho;
    jac.set(Co55, H1, scratch);

    scratch = -screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(Co55)*state.rho - screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*state.rho + screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(Fe52)*state.rho + screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(Mn51)*state.rho + screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(Mn52)*state.rho;
    jac.set(Co55, He4, scratch);

    jac.set(Co55, C12, 0.0);

    jac.set(Co55, C13, 0.0);

    jac.set(Co55, C14, 0.0);

    jac.set(Co55, N13, 0.0);

    jac.set(Co55, N14, 0.0);

    jac.set(Co55, N15, 0.0);

    jac.set(Co55, O16, 0.0);

    jac.set(Co55, O17, 0.0);

    jac.set(Co55, O18, 0.0);

    jac.set(Co55, F18, 0.0);

    jac.set(Co55, F19, 0.0);

    jac.set(Co55, Ne20, 0.0);

    jac.set(Co55, Ne21, 0.0);

    jac.set(Co55, Ne22, 0.0);

    jac.set(Co55, Na22, 0.0);

    jac.set(Co55, Na23, 0.0);

    jac.set(Co55, Mg24, 0.0);

    jac.set(Co55, Mg25, 0.0);

    jac.set(Co55, Al26, 0.0);

    jac.set(Co55, Al27, 0.0);

    jac.set(Co55, Si28, 0.0);

    jac.set(Co55, Si29, 0.0);

    jac.set(Co55, P30, 0.0);

    jac.set(Co55, P31, 0.0);

    jac.set(Co55, S32, 0.0);

    jac.set(Co55, S33, 0.0);

    jac.set(Co55, Cl34, 0.0);

    jac.set(Co55, Cl35, 0.0);

    jac.set(Co55, Ar36, 0.0);

    jac.set(Co55, Ar37, 0.0);

    jac.set(Co55, K38, 0.0);

    jac.set(Co55, K39, 0.0);

    jac.set(Co55, Ca40, 0.0);

    jac.set(Co55, Ca41, 0.0);

    jac.set(Co55, Sc42, 0.0);

    jac.set(Co55, Sc43, 0.0);

    jac.set(Co55, Ti44, 0.0);

    jac.set(Co55, Ti45, 0.0);

    jac.set(Co55, V46, 0.0);

    jac.set(Co55, V47, 0.0);

    jac.set(Co55, Cr48, 0.0);

    jac.set(Co55, Cr49, 0.0);

    jac.set(Co55, Mn50, 0.0);

    scratch = screened_rates(k_He4_Mn51_to_Co55_reaclib)*Y(He4)*state.rho;
    jac.set(Co55, Mn51, scratch);

    scratch = screened_rates(k_He4_Mn52_to_n_Co55_derived)*Y(He4)*state.rho;
    jac.set(Co55, Mn52, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_Co55_reaclib)*Y(He4)*state.rho;
    jac.set(Co55, Fe52, scratch);

    jac.set(Co55, Fe53, 0.0);

    scratch = screened_rates(k_p_Fe54_to_Co55_reaclib)*Y(H1)*state.rho;
    jac.set(Co55, Fe54, scratch);

    scratch = screened_rates(k_Fe55_to_Co55_weaktab) + screened_rates(k_p_Fe55_to_n_Co55_derived)*Y(H1)*state.rho;
    jac.set(Co55, Fe55, scratch);

    jac.set(Co55, Fe56, 0.0);

    scratch = -screened_rates(k_Co55_to_Fe55_weaktab) - screened_rates(k_Co55_to_He4_Mn51_derived) - screened_rates(k_Co55_to_p_Fe54_derived) - screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(He4)*state.rho - screened_rates(k_n_Co55_to_Co56_reaclib)*Y(N)*state.rho - screened_rates(k_n_Co55_to_He4_Mn52_reaclib)*Y(N)*state.rho - screened_rates(k_n_Co55_to_p_Fe55_reaclib)*Y(N)*state.rho - screened_rates(k_p_Co55_to_He4_Fe52_derived)*Y(H1)*state.rho - screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(Co55, Co55, scratch);

    scratch = screened_rates(k_Co56_to_n_Co55_derived);
    jac.set(Co55, Co56, scratch);

    jac.set(Co55, Co57, 0.0);

    scratch = screened_rates(k_Ni56_to_p_Co55_derived);
    jac.set(Co55, Ni56, scratch);

    jac.set(Co55, Ni57, 0.0);

    scratch = screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(H1)*state.rho;
    jac.set(Co55, Ni58, scratch);

    scratch = screened_rates(k_Cu59_to_He4_Co55_derived);
    jac.set(Co55, Cu59, scratch);

    jac.set(Co55, Zn60, 0.0);

    scratch = screened_rates(k_n_Co55_to_Co56_reaclib)*Y(Co55)*state.rho - screened_rates(k_n_Co56_to_Co57_reaclib)*Y(Co56)*state.rho - screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(Co56)*state.rho + screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*state.rho + screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*state.rho;
    jac.set(Co56, N, scratch);

    scratch = -screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(Co56)*state.rho - screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(Co56)*state.rho - screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho + screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(Fe55)*state.rho + screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho;
    jac.set(Co56, H1, scratch);

    scratch = -screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*state.rho + screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(Fe53)*state.rho + screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(Mn52)*state.rho;
    jac.set(Co56, He4, scratch);

    jac.set(Co56, C12, 0.0);

    jac.set(Co56, C13, 0.0);

    jac.set(Co56, C14, 0.0);

    jac.set(Co56, N13, 0.0);

    jac.set(Co56, N14, 0.0);

    jac.set(Co56, N15, 0.0);

    jac.set(Co56, O16, 0.0);

    jac.set(Co56, O17, 0.0);

    jac.set(Co56, O18, 0.0);

    jac.set(Co56, F18, 0.0);

    jac.set(Co56, F19, 0.0);

    jac.set(Co56, Ne20, 0.0);

    jac.set(Co56, Ne21, 0.0);

    jac.set(Co56, Ne22, 0.0);

    jac.set(Co56, Na22, 0.0);

    jac.set(Co56, Na23, 0.0);

    jac.set(Co56, Mg24, 0.0);

    jac.set(Co56, Mg25, 0.0);

    jac.set(Co56, Al26, 0.0);

    jac.set(Co56, Al27, 0.0);

    jac.set(Co56, Si28, 0.0);

    jac.set(Co56, Si29, 0.0);

    jac.set(Co56, P30, 0.0);

    jac.set(Co56, P31, 0.0);

    jac.set(Co56, S32, 0.0);

    jac.set(Co56, S33, 0.0);

    jac.set(Co56, Cl34, 0.0);

    jac.set(Co56, Cl35, 0.0);

    jac.set(Co56, Ar36, 0.0);

    jac.set(Co56, Ar37, 0.0);

    jac.set(Co56, K38, 0.0);

    jac.set(Co56, K39, 0.0);

    jac.set(Co56, Ca40, 0.0);

    jac.set(Co56, Ca41, 0.0);

    jac.set(Co56, Sc42, 0.0);

    jac.set(Co56, Sc43, 0.0);

    jac.set(Co56, Ti44, 0.0);

    jac.set(Co56, Ti45, 0.0);

    jac.set(Co56, V46, 0.0);

    jac.set(Co56, V47, 0.0);

    jac.set(Co56, Cr48, 0.0);

    jac.set(Co56, Cr49, 0.0);

    jac.set(Co56, Mn50, 0.0);

    jac.set(Co56, Mn51, 0.0);

    scratch = screened_rates(k_He4_Mn52_to_Co56_reaclib)*Y(He4)*state.rho;
    jac.set(Co56, Mn52, scratch);

    jac.set(Co56, Fe52, 0.0);

    scratch = screened_rates(k_He4_Fe53_to_p_Co56_reaclib)*Y(He4)*state.rho;
    jac.set(Co56, Fe53, scratch);

    jac.set(Co56, Fe54, 0.0);

    scratch = screened_rates(k_p_Fe55_to_Co56_reaclib)*Y(H1)*state.rho;
    jac.set(Co56, Fe55, scratch);

    scratch = screened_rates(k_Fe56_to_Co56_weaktab) + screened_rates(k_p_Fe56_to_n_Co56_derived)*Y(H1)*state.rho;
    jac.set(Co56, Fe56, scratch);

    scratch = screened_rates(k_n_Co55_to_Co56_reaclib)*Y(N)*state.rho;
    jac.set(Co56, Co55, scratch);

    scratch = -screened_rates(k_Co56_to_Fe56_weaktab) - screened_rates(k_Co56_to_He4_Mn52_derived) - screened_rates(k_Co56_to_Ni56_weaktab) - screened_rates(k_Co56_to_n_Co55_derived) - screened_rates(k_Co56_to_p_Fe55_derived) - screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(He4)*state.rho - screened_rates(k_n_Co56_to_Co57_reaclib)*Y(N)*state.rho - screened_rates(k_n_Co56_to_p_Fe56_reaclib)*Y(N)*state.rho - screened_rates(k_p_Co56_to_He4_Fe53_derived)*Y(H1)*state.rho - screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(H1)*state.rho;
    jac.set(Co56, Co56, scratch);

    scratch = screened_rates(k_Co57_to_n_Co56_derived);
    jac.set(Co56, Co57, scratch);

    scratch = screened_rates(k_Ni56_to_Co56_weaktab) + screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(N)*state.rho;
    jac.set(Co56, Ni56, scratch);

    scratch = screened_rates(k_Ni57_to_p_Co56_derived);
    jac.set(Co56, Ni57, scratch);

    jac.set(Co56, Ni58, 0.0);

    scratch = screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(N)*state.rho;
    jac.set(Co56, Cu59, scratch);

    jac.set(Co56, Zn60, 0.0);

    scratch = screened_rates(k_n_Co56_to_Co57_reaclib)*Y(Co56)*state.rho + screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*state.rho;
    jac.set(Co57, N, scratch);

    scratch = -screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(Co57)*state.rho - screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(Co57)*state.rho - screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*state.rho + screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(Fe56)*state.rho;
    jac.set(Co57, H1, scratch);

    scratch = screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(Fe54)*state.rho;
    jac.set(Co57, He4, scratch);

    jac.set(Co57, C12, 0.0);

    jac.set(Co57, C13, 0.0);

    jac.set(Co57, C14, 0.0);

    jac.set(Co57, N13, 0.0);

    jac.set(Co57, N14, 0.0);

    jac.set(Co57, N15, 0.0);

    jac.set(Co57, O16, 0.0);

    jac.set(Co57, O17, 0.0);

    jac.set(Co57, O18, 0.0);

    jac.set(Co57, F18, 0.0);

    jac.set(Co57, F19, 0.0);

    jac.set(Co57, Ne20, 0.0);

    jac.set(Co57, Ne21, 0.0);

    jac.set(Co57, Ne22, 0.0);

    jac.set(Co57, Na22, 0.0);

    jac.set(Co57, Na23, 0.0);

    jac.set(Co57, Mg24, 0.0);

    jac.set(Co57, Mg25, 0.0);

    jac.set(Co57, Al26, 0.0);

    jac.set(Co57, Al27, 0.0);

    jac.set(Co57, Si28, 0.0);

    jac.set(Co57, Si29, 0.0);

    jac.set(Co57, P30, 0.0);

    jac.set(Co57, P31, 0.0);

    jac.set(Co57, S32, 0.0);

    jac.set(Co57, S33, 0.0);

    jac.set(Co57, Cl34, 0.0);

    jac.set(Co57, Cl35, 0.0);

    jac.set(Co57, Ar36, 0.0);

    jac.set(Co57, Ar37, 0.0);

    jac.set(Co57, K38, 0.0);

    jac.set(Co57, K39, 0.0);

    jac.set(Co57, Ca40, 0.0);

    jac.set(Co57, Ca41, 0.0);

    jac.set(Co57, Sc42, 0.0);

    jac.set(Co57, Sc43, 0.0);

    jac.set(Co57, Ti44, 0.0);

    jac.set(Co57, Ti45, 0.0);

    jac.set(Co57, V46, 0.0);

    jac.set(Co57, V47, 0.0);

    jac.set(Co57, Cr48, 0.0);

    jac.set(Co57, Cr49, 0.0);

    jac.set(Co57, Mn50, 0.0);

    jac.set(Co57, Mn51, 0.0);

    jac.set(Co57, Mn52, 0.0);

    jac.set(Co57, Fe52, 0.0);

    jac.set(Co57, Fe53, 0.0);

    scratch = screened_rates(k_He4_Fe54_to_p_Co57_derived)*Y(He4)*state.rho;
    jac.set(Co57, Fe54, scratch);

    jac.set(Co57, Fe55, 0.0);

    scratch = screened_rates(k_p_Fe56_to_Co57_reaclib)*Y(H1)*state.rho;
    jac.set(Co57, Fe56, scratch);

    jac.set(Co57, Co55, 0.0);

    scratch = screened_rates(k_n_Co56_to_Co57_reaclib)*Y(N)*state.rho;
    jac.set(Co57, Co56, scratch);

    scratch = -screened_rates(k_Co57_to_Ni57_weaktab) - screened_rates(k_Co57_to_n_Co56_derived) - screened_rates(k_Co57_to_p_Fe56_derived) - screened_rates(k_p_Co57_to_He4_Fe54_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(H1)*state.rho;
    jac.set(Co57, Co57, scratch);

    jac.set(Co57, Ni56, 0.0);

    scratch = screened_rates(k_Ni57_to_Co57_weaktab) + screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(N)*state.rho;
    jac.set(Co57, Ni57, scratch);

    scratch = screened_rates(k_Ni58_to_p_Co57_derived);
    jac.set(Co57, Ni58, scratch);

    jac.set(Co57, Cu59, 0.0);

    jac.set(Co57, Zn60, 0.0);

    scratch = -screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(Ni56)*state.rho - screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(Ni56)*state.rho - screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(Ni56)*state.rho;
    jac.set(Ni56, N, scratch);

    scratch = screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(Co55)*state.rho + screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho + screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho;
    jac.set(Ni56, H1, scratch);

    scratch = screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(Fe52)*state.rho + screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(Fe53)*state.rho - screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(Ni56)*state.rho - screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(Ni56)*state.rho;
    jac.set(Ni56, He4, scratch);

    jac.set(Ni56, C12, 0.0);

    jac.set(Ni56, C13, 0.0);

    jac.set(Ni56, C14, 0.0);

    jac.set(Ni56, N13, 0.0);

    jac.set(Ni56, N14, 0.0);

    jac.set(Ni56, N15, 0.0);

    jac.set(Ni56, O16, 0.0);

    jac.set(Ni56, O17, 0.0);

    jac.set(Ni56, O18, 0.0);

    jac.set(Ni56, F18, 0.0);

    jac.set(Ni56, F19, 0.0);

    jac.set(Ni56, Ne20, 0.0);

    jac.set(Ni56, Ne21, 0.0);

    jac.set(Ni56, Ne22, 0.0);

    jac.set(Ni56, Na22, 0.0);

    jac.set(Ni56, Na23, 0.0);

    jac.set(Ni56, Mg24, 0.0);

    jac.set(Ni56, Mg25, 0.0);

    jac.set(Ni56, Al26, 0.0);

    jac.set(Ni56, Al27, 0.0);

    jac.set(Ni56, Si28, 0.0);

    jac.set(Ni56, Si29, 0.0);

    jac.set(Ni56, P30, 0.0);

    jac.set(Ni56, P31, 0.0);

    jac.set(Ni56, S32, 0.0);

    jac.set(Ni56, S33, 0.0);

    jac.set(Ni56, Cl34, 0.0);

    jac.set(Ni56, Cl35, 0.0);

    jac.set(Ni56, Ar36, 0.0);

    jac.set(Ni56, Ar37, 0.0);

    jac.set(Ni56, K38, 0.0);

    jac.set(Ni56, K39, 0.0);

    jac.set(Ni56, Ca40, 0.0);

    jac.set(Ni56, Ca41, 0.0);

    jac.set(Ni56, Sc42, 0.0);

    jac.set(Ni56, Sc43, 0.0);

    jac.set(Ni56, Ti44, 0.0);

    jac.set(Ni56, Ti45, 0.0);

    jac.set(Ni56, V46, 0.0);

    jac.set(Ni56, V47, 0.0);

    jac.set(Ni56, Cr48, 0.0);

    jac.set(Ni56, Cr49, 0.0);

    jac.set(Ni56, Mn50, 0.0);

    jac.set(Ni56, Mn51, 0.0);

    jac.set(Ni56, Mn52, 0.0);

    scratch = screened_rates(k_He4_Fe52_to_Ni56_reaclib)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = screened_rates(k_He4_Fe53_to_n_Ni56_derived)*Y(He4)*state.rho;
    jac.set(Ni56, Fe53, scratch);

    jac.set(Ni56, Fe54, 0.0);

    jac.set(Ni56, Fe55, 0.0);

    jac.set(Ni56, Fe56, 0.0);

    scratch = screened_rates(k_p_Co55_to_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(Ni56, Co55, scratch);

    scratch = screened_rates(k_Co56_to_Ni56_weaktab) + screened_rates(k_p_Co56_to_n_Ni56_derived)*Y(H1)*state.rho;
    jac.set(Ni56, Co56, scratch);

    jac.set(Ni56, Co57, 0.0);

    scratch = -screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*state.rho - screened_rates(k_Ni56_to_Co56_weaktab) - screened_rates(k_Ni56_to_He4_Fe52_derived) - screened_rates(k_Ni56_to_p_Co55_derived) - screened_rates(k_n_Ni56_to_He4_Fe53_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni56_to_p_Co56_reaclib)*Y(N)*state.rho;
    jac.set(Ni56, Ni56, scratch);

    scratch = screened_rates(k_Ni57_to_n_Ni56_derived);
    jac.set(Ni56, Ni57, scratch);

    jac.set(Ni56, Ni58, 0.0);

    scratch = screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(Ni56, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_He4_Ni56_derived);
    jac.set(Ni56, Zn60, scratch);

    scratch = screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(Ni56)*state.rho - screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(Ni57)*state.rho - screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(Ni57)*state.rho - screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(Ni57)*state.rho + screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*state.rho;
    jac.set(Ni57, N, scratch);

    scratch = screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(Co56)*state.rho + screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(Co57)*state.rho;
    jac.set(Ni57, H1, scratch);

    scratch = screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(Fe53)*state.rho + screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(Fe54)*state.rho - screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(Ni57)*state.rho;
    jac.set(Ni57, He4, scratch);

    jac.set(Ni57, C12, 0.0);

    jac.set(Ni57, C13, 0.0);

    jac.set(Ni57, C14, 0.0);

    jac.set(Ni57, N13, 0.0);

    jac.set(Ni57, N14, 0.0);

    jac.set(Ni57, N15, 0.0);

    jac.set(Ni57, O16, 0.0);

    jac.set(Ni57, O17, 0.0);

    jac.set(Ni57, O18, 0.0);

    jac.set(Ni57, F18, 0.0);

    jac.set(Ni57, F19, 0.0);

    jac.set(Ni57, Ne20, 0.0);

    jac.set(Ni57, Ne21, 0.0);

    jac.set(Ni57, Ne22, 0.0);

    jac.set(Ni57, Na22, 0.0);

    jac.set(Ni57, Na23, 0.0);

    jac.set(Ni57, Mg24, 0.0);

    jac.set(Ni57, Mg25, 0.0);

    jac.set(Ni57, Al26, 0.0);

    jac.set(Ni57, Al27, 0.0);

    jac.set(Ni57, Si28, 0.0);

    jac.set(Ni57, Si29, 0.0);

    jac.set(Ni57, P30, 0.0);

    jac.set(Ni57, P31, 0.0);

    jac.set(Ni57, S32, 0.0);

    jac.set(Ni57, S33, 0.0);

    jac.set(Ni57, Cl34, 0.0);

    jac.set(Ni57, Cl35, 0.0);

    jac.set(Ni57, Ar36, 0.0);

    jac.set(Ni57, Ar37, 0.0);

    jac.set(Ni57, K38, 0.0);

    jac.set(Ni57, K39, 0.0);

    jac.set(Ni57, Ca40, 0.0);

    jac.set(Ni57, Ca41, 0.0);

    jac.set(Ni57, Sc42, 0.0);

    jac.set(Ni57, Sc43, 0.0);

    jac.set(Ni57, Ti44, 0.0);

    jac.set(Ni57, Ti45, 0.0);

    jac.set(Ni57, V46, 0.0);

    jac.set(Ni57, V47, 0.0);

    jac.set(Ni57, Cr48, 0.0);

    jac.set(Ni57, Cr49, 0.0);

    jac.set(Ni57, Mn50, 0.0);

    jac.set(Ni57, Mn51, 0.0);

    jac.set(Ni57, Mn52, 0.0);

    jac.set(Ni57, Fe52, 0.0);

    scratch = screened_rates(k_He4_Fe53_to_Ni57_reaclib)*Y(He4)*state.rho;
    jac.set(Ni57, Fe53, scratch);

    scratch = screened_rates(k_He4_Fe54_to_n_Ni57_derived)*Y(He4)*state.rho;
    jac.set(Ni57, Fe54, scratch);

    jac.set(Ni57, Fe55, 0.0);

    jac.set(Ni57, Fe56, 0.0);

    jac.set(Ni57, Co55, 0.0);

    scratch = screened_rates(k_p_Co56_to_Ni57_reaclib)*Y(H1)*state.rho;
    jac.set(Ni57, Co56, scratch);

    scratch = screened_rates(k_Co57_to_Ni57_weaktab) + screened_rates(k_p_Co57_to_n_Ni57_derived)*Y(H1)*state.rho;
    jac.set(Ni57, Co57, scratch);

    scratch = screened_rates(k_n_Ni56_to_Ni57_reaclib)*Y(N)*state.rho;
    jac.set(Ni57, Ni56, scratch);

    scratch = -screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*state.rho - screened_rates(k_Ni57_to_Co57_weaktab) - screened_rates(k_Ni57_to_He4_Fe53_derived) - screened_rates(k_Ni57_to_n_Ni56_derived) - screened_rates(k_Ni57_to_p_Co56_derived) - screened_rates(k_n_Ni57_to_He4_Fe54_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(N)*state.rho - screened_rates(k_n_Ni57_to_p_Co57_reaclib)*Y(N)*state.rho;
    jac.set(Ni57, Ni57, scratch);

    scratch = screened_rates(k_Ni58_to_n_Ni57_derived);
    jac.set(Ni57, Ni58, scratch);

    jac.set(Ni57, Cu59, 0.0);

    scratch = screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(N)*state.rho;
    jac.set(Ni57, Zn60, scratch);

    scratch = screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(Ni57)*state.rho - screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(Ni58)*state.rho;
    jac.set(Ni58, N, scratch);

    scratch = screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(Co57)*state.rho - screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(Ni58)*state.rho - screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho;
    jac.set(Ni58, H1, scratch);

    scratch = screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(Co55)*state.rho + screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(Fe54)*state.rho + screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(Fe55)*state.rho;
    jac.set(Ni58, He4, scratch);

    jac.set(Ni58, C12, 0.0);

    jac.set(Ni58, C13, 0.0);

    jac.set(Ni58, C14, 0.0);

    jac.set(Ni58, N13, 0.0);

    jac.set(Ni58, N14, 0.0);

    jac.set(Ni58, N15, 0.0);

    jac.set(Ni58, O16, 0.0);

    jac.set(Ni58, O17, 0.0);

    jac.set(Ni58, O18, 0.0);

    jac.set(Ni58, F18, 0.0);

    jac.set(Ni58, F19, 0.0);

    jac.set(Ni58, Ne20, 0.0);

    jac.set(Ni58, Ne21, 0.0);

    jac.set(Ni58, Ne22, 0.0);

    jac.set(Ni58, Na22, 0.0);

    jac.set(Ni58, Na23, 0.0);

    jac.set(Ni58, Mg24, 0.0);

    jac.set(Ni58, Mg25, 0.0);

    jac.set(Ni58, Al26, 0.0);

    jac.set(Ni58, Al27, 0.0);

    jac.set(Ni58, Si28, 0.0);

    jac.set(Ni58, Si29, 0.0);

    jac.set(Ni58, P30, 0.0);

    jac.set(Ni58, P31, 0.0);

    jac.set(Ni58, S32, 0.0);

    jac.set(Ni58, S33, 0.0);

    jac.set(Ni58, Cl34, 0.0);

    jac.set(Ni58, Cl35, 0.0);

    jac.set(Ni58, Ar36, 0.0);

    jac.set(Ni58, Ar37, 0.0);

    jac.set(Ni58, K38, 0.0);

    jac.set(Ni58, K39, 0.0);

    jac.set(Ni58, Ca40, 0.0);

    jac.set(Ni58, Ca41, 0.0);

    jac.set(Ni58, Sc42, 0.0);

    jac.set(Ni58, Sc43, 0.0);

    jac.set(Ni58, Ti44, 0.0);

    jac.set(Ni58, Ti45, 0.0);

    jac.set(Ni58, V46, 0.0);

    jac.set(Ni58, V47, 0.0);

    jac.set(Ni58, Cr48, 0.0);

    jac.set(Ni58, Cr49, 0.0);

    jac.set(Ni58, Mn50, 0.0);

    jac.set(Ni58, Mn51, 0.0);

    jac.set(Ni58, Mn52, 0.0);

    jac.set(Ni58, Fe52, 0.0);

    jac.set(Ni58, Fe53, 0.0);

    scratch = screened_rates(k_He4_Fe54_to_Ni58_reaclib)*Y(He4)*state.rho;
    jac.set(Ni58, Fe54, scratch);

    scratch = screened_rates(k_He4_Fe55_to_n_Ni58_derived)*Y(He4)*state.rho;
    jac.set(Ni58, Fe55, scratch);

    jac.set(Ni58, Fe56, 0.0);

    scratch = screened_rates(k_He4_Co55_to_p_Ni58_reaclib)*Y(He4)*state.rho;
    jac.set(Ni58, Co55, scratch);

    jac.set(Ni58, Co56, 0.0);

    scratch = screened_rates(k_p_Co57_to_Ni58_reaclib)*Y(H1)*state.rho;
    jac.set(Ni58, Co57, scratch);

    jac.set(Ni58, Ni56, 0.0);

    scratch = screened_rates(k_n_Ni57_to_Ni58_reaclib)*Y(N)*state.rho;
    jac.set(Ni58, Ni57, scratch);

    scratch = -screened_rates(k_Ni58_to_He4_Fe54_derived) - screened_rates(k_Ni58_to_n_Ni57_derived) - screened_rates(k_Ni58_to_p_Co57_derived) - screened_rates(k_n_Ni58_to_He4_Fe55_reaclib)*Y(N)*state.rho - screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Ni58_to_He4_Co55_derived)*Y(H1)*state.rho;
    jac.set(Ni58, Ni58, scratch);

    scratch = screened_rates(k_Cu59_to_p_Ni58_derived);
    jac.set(Ni58, Cu59, scratch);

    jac.set(Ni58, Zn60, 0.0);

    scratch = -screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(Cu59)*state.rho;
    jac.set(Cu59, N, scratch);

    scratch = -screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho - screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(Cu59)*state.rho + screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(Ni58)*state.rho;
    jac.set(Cu59, H1, scratch);

    scratch = screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(Co55)*state.rho + screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(Co56)*state.rho + screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(Ni56)*state.rho;
    jac.set(Cu59, He4, scratch);

    jac.set(Cu59, C12, 0.0);

    jac.set(Cu59, C13, 0.0);

    jac.set(Cu59, C14, 0.0);

    jac.set(Cu59, N13, 0.0);

    jac.set(Cu59, N14, 0.0);

    jac.set(Cu59, N15, 0.0);

    jac.set(Cu59, O16, 0.0);

    jac.set(Cu59, O17, 0.0);

    jac.set(Cu59, O18, 0.0);

    jac.set(Cu59, F18, 0.0);

    jac.set(Cu59, F19, 0.0);

    jac.set(Cu59, Ne20, 0.0);

    jac.set(Cu59, Ne21, 0.0);

    jac.set(Cu59, Ne22, 0.0);

    jac.set(Cu59, Na22, 0.0);

    jac.set(Cu59, Na23, 0.0);

    jac.set(Cu59, Mg24, 0.0);

    jac.set(Cu59, Mg25, 0.0);

    jac.set(Cu59, Al26, 0.0);

    jac.set(Cu59, Al27, 0.0);

    jac.set(Cu59, Si28, 0.0);

    jac.set(Cu59, Si29, 0.0);

    jac.set(Cu59, P30, 0.0);

    jac.set(Cu59, P31, 0.0);

    jac.set(Cu59, S32, 0.0);

    jac.set(Cu59, S33, 0.0);

    jac.set(Cu59, Cl34, 0.0);

    jac.set(Cu59, Cl35, 0.0);

    jac.set(Cu59, Ar36, 0.0);

    jac.set(Cu59, Ar37, 0.0);

    jac.set(Cu59, K38, 0.0);

    jac.set(Cu59, K39, 0.0);

    jac.set(Cu59, Ca40, 0.0);

    jac.set(Cu59, Ca41, 0.0);

    jac.set(Cu59, Sc42, 0.0);

    jac.set(Cu59, Sc43, 0.0);

    jac.set(Cu59, Ti44, 0.0);

    jac.set(Cu59, Ti45, 0.0);

    jac.set(Cu59, V46, 0.0);

    jac.set(Cu59, V47, 0.0);

    jac.set(Cu59, Cr48, 0.0);

    jac.set(Cu59, Cr49, 0.0);

    jac.set(Cu59, Mn50, 0.0);

    jac.set(Cu59, Mn51, 0.0);

    jac.set(Cu59, Mn52, 0.0);

    jac.set(Cu59, Fe52, 0.0);

    jac.set(Cu59, Fe53, 0.0);

    jac.set(Cu59, Fe54, 0.0);

    jac.set(Cu59, Fe55, 0.0);

    jac.set(Cu59, Fe56, 0.0);

    scratch = screened_rates(k_He4_Co55_to_Cu59_reaclib)*Y(He4)*state.rho;
    jac.set(Cu59, Co55, scratch);

    scratch = screened_rates(k_He4_Co56_to_n_Cu59_derived)*Y(He4)*state.rho;
    jac.set(Cu59, Co56, scratch);

    jac.set(Cu59, Co57, 0.0);

    scratch = screened_rates(k_He4_Ni56_to_p_Cu59_derived)*Y(He4)*state.rho;
    jac.set(Cu59, Ni56, scratch);

    jac.set(Cu59, Ni57, 0.0);

    scratch = screened_rates(k_p_Ni58_to_Cu59_reaclib)*Y(H1)*state.rho;
    jac.set(Cu59, Ni58, scratch);

    scratch = -screened_rates(k_Cu59_to_He4_Co55_derived) - screened_rates(k_Cu59_to_p_Ni58_derived) - screened_rates(k_n_Cu59_to_He4_Co56_reaclib)*Y(N)*state.rho - screened_rates(k_p_Cu59_to_He4_Ni56_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(H1)*state.rho;
    jac.set(Cu59, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_p_Cu59_derived);
    jac.set(Cu59, Zn60, scratch);

    scratch = -screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(Zn60)*state.rho;
    jac.set(Zn60, N, scratch);

    scratch = screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(Cu59)*state.rho;
    jac.set(Zn60, H1, scratch);

    scratch = screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(Ni56)*state.rho + screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(Ni57)*state.rho;
    jac.set(Zn60, He4, scratch);

    jac.set(Zn60, C12, 0.0);

    jac.set(Zn60, C13, 0.0);

    jac.set(Zn60, C14, 0.0);

    jac.set(Zn60, N13, 0.0);

    jac.set(Zn60, N14, 0.0);

    jac.set(Zn60, N15, 0.0);

    jac.set(Zn60, O16, 0.0);

    jac.set(Zn60, O17, 0.0);

    jac.set(Zn60, O18, 0.0);

    jac.set(Zn60, F18, 0.0);

    jac.set(Zn60, F19, 0.0);

    jac.set(Zn60, Ne20, 0.0);

    jac.set(Zn60, Ne21, 0.0);

    jac.set(Zn60, Ne22, 0.0);

    jac.set(Zn60, Na22, 0.0);

    jac.set(Zn60, Na23, 0.0);

    jac.set(Zn60, Mg24, 0.0);

    jac.set(Zn60, Mg25, 0.0);

    jac.set(Zn60, Al26, 0.0);

    jac.set(Zn60, Al27, 0.0);

    jac.set(Zn60, Si28, 0.0);

    jac.set(Zn60, Si29, 0.0);

    jac.set(Zn60, P30, 0.0);

    jac.set(Zn60, P31, 0.0);

    jac.set(Zn60, S32, 0.0);

    jac.set(Zn60, S33, 0.0);

    jac.set(Zn60, Cl34, 0.0);

    jac.set(Zn60, Cl35, 0.0);

    jac.set(Zn60, Ar36, 0.0);

    jac.set(Zn60, Ar37, 0.0);

    jac.set(Zn60, K38, 0.0);

    jac.set(Zn60, K39, 0.0);

    jac.set(Zn60, Ca40, 0.0);

    jac.set(Zn60, Ca41, 0.0);

    jac.set(Zn60, Sc42, 0.0);

    jac.set(Zn60, Sc43, 0.0);

    jac.set(Zn60, Ti44, 0.0);

    jac.set(Zn60, Ti45, 0.0);

    jac.set(Zn60, V46, 0.0);

    jac.set(Zn60, V47, 0.0);

    jac.set(Zn60, Cr48, 0.0);

    jac.set(Zn60, Cr49, 0.0);

    jac.set(Zn60, Mn50, 0.0);

    jac.set(Zn60, Mn51, 0.0);

    jac.set(Zn60, Mn52, 0.0);

    jac.set(Zn60, Fe52, 0.0);

    jac.set(Zn60, Fe53, 0.0);

    jac.set(Zn60, Fe54, 0.0);

    jac.set(Zn60, Fe55, 0.0);

    jac.set(Zn60, Fe56, 0.0);

    jac.set(Zn60, Co55, 0.0);

    jac.set(Zn60, Co56, 0.0);

    jac.set(Zn60, Co57, 0.0);

    scratch = screened_rates(k_He4_Ni56_to_Zn60_reaclib)*Y(He4)*state.rho;
    jac.set(Zn60, Ni56, scratch);

    scratch = screened_rates(k_He4_Ni57_to_n_Zn60_derived)*Y(He4)*state.rho;
    jac.set(Zn60, Ni57, scratch);

    jac.set(Zn60, Ni58, 0.0);

    scratch = screened_rates(k_p_Cu59_to_Zn60_reaclib)*Y(H1)*state.rho;
    jac.set(Zn60, Cu59, scratch);

    scratch = -screened_rates(k_Zn60_to_He4_Ni56_derived) - screened_rates(k_Zn60_to_p_Cu59_derived) - screened_rates(k_n_Zn60_to_He4_Ni57_reaclib)*Y(N)*state.rho;
    jac.set(Zn60, Zn60, scratch);


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

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> amrex::Real { return jac.get(i, j); };
        amrex::Real tmp{};
        ener_gener_rate(jac_slice_2, tmp);
        jac(net_ienuc,j) = tmp;
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

    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    amrex::Array1D<amrex::Real, 1, neqs>  yderivs;

    rhs_nuc(state, yderivs, Y, rate_eval.dscreened_rates_dT);

    for (int k = 1; k <= NumSpec; k++) {
        jac.set(k, net_ienuc, temperature_to_energy_jacobian(state, yderivs(k)));
    }


    // finally, d(de/dt)/de

    amrex::Real jac_e_T;
    ener_gener_rate(yderivs, jac_e_T);
    jac_e_T -= dsneutdt;
    jac.set(net_ienuc, net_ienuc, temperature_to_energy_jacobian(state, jac_e_T));

}


AMREX_INLINE
void actual_rhs_init () {


}


#endif
