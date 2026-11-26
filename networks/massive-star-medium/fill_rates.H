#ifndef FILL_RATES_H
#define FILL_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

namespace rate_util {

    using namespace Rates;
    using namespace Species;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] const amrex::Real temp,
               [[maybe_unused]] const amrex::Real rho,
               [[maybe_unused]] const amrex::Real rhoy,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate{};
        [[maybe_unused]] amrex::Real drate_dT{};

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

        // rates are filled in a topologically-sorted order

        {
            // O16_to_He4_C12_reaclib
            reaclib_rates::rate_O16_to_He4_C12_reaclib(rate_eval, tfactors);
        }

        {
            // Ne20_to_He4_O16_reaclib
            reaclib_rates::rate_Ne20_to_He4_O16_reaclib(rate_eval, tfactors);
        }

        {
            // Mg24_to_n_Mg23_reaclib
            reaclib_rates::rate_Mg24_to_n_Mg23_reaclib(rate_eval, tfactors);
        }

        {
            // Si27_to_He4_Mg23_reaclib
            reaclib_rates::rate_Si27_to_He4_Mg23_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_to_n_Si27_reaclib
            reaclib_rates::rate_Si28_to_n_Si27_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_to_p_Al27_reaclib
            reaclib_rates::rate_Si28_to_p_Al27_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_to_He4_Mg24_reaclib
            reaclib_rates::rate_Si28_to_He4_Mg24_reaclib(rate_eval, tfactors);
        }

        {
            // Si29_to_n_Si28_reaclib
            reaclib_rates::rate_Si29_to_n_Si28_reaclib(rate_eval, tfactors);
        }

        {
            // Si30_to_n_Si29_reaclib
            reaclib_rates::rate_Si30_to_n_Si29_reaclib(rate_eval, tfactors);
        }

        {
            // Si31_to_n_Si30_reaclib
            reaclib_rates::rate_Si31_to_n_Si30_reaclib(rate_eval, tfactors);
        }

        {
            // P29_to_p_Si28_reaclib
            reaclib_rates::rate_P29_to_p_Si28_reaclib(rate_eval, tfactors);
        }

        {
            // P31_to_p_Si30_reaclib
            reaclib_rates::rate_P31_to_p_Si30_reaclib(rate_eval, tfactors);
        }

        {
            // P31_to_He4_Al27_reaclib
            reaclib_rates::rate_P31_to_He4_Al27_reaclib(rate_eval, tfactors);
        }

        {
            // S32_to_p_P31_reaclib
            reaclib_rates::rate_S32_to_p_P31_reaclib(rate_eval, tfactors);
        }

        {
            // S32_to_He4_Si28_reaclib
            reaclib_rates::rate_S32_to_He4_Si28_reaclib(rate_eval, tfactors);
        }

        {
            // S33_to_n_S32_reaclib
            reaclib_rates::rate_S33_to_n_S32_reaclib(rate_eval, tfactors);
        }

        {
            // S33_to_He4_Si29_reaclib
            reaclib_rates::rate_S33_to_He4_Si29_reaclib(rate_eval, tfactors);
        }

        {
            // S34_to_n_S33_reaclib
            reaclib_rates::rate_S34_to_n_S33_reaclib(rate_eval, tfactors);
        }

        {
            // S34_to_p_P33_reaclib
            reaclib_rates::rate_S34_to_p_P33_reaclib(rate_eval, tfactors);
        }

        {
            // S34_to_He4_Si30_reaclib
            reaclib_rates::rate_S34_to_He4_Si30_reaclib(rate_eval, tfactors);
        }

        {
            // S35_to_n_S34_reaclib
            reaclib_rates::rate_S35_to_n_S34_reaclib(rate_eval, tfactors);
        }

        {
            // S35_to_He4_Si31_reaclib
            reaclib_rates::rate_S35_to_He4_Si31_reaclib(rate_eval, tfactors);
        }

        {
            // Cl33_to_p_S32_reaclib
            reaclib_rates::rate_Cl33_to_p_S32_reaclib(rate_eval, tfactors);
        }

        {
            // Cl33_to_He4_P29_reaclib
            reaclib_rates::rate_Cl33_to_He4_P29_reaclib(rate_eval, tfactors);
        }

        {
            // Cl35_to_p_S34_reaclib
            reaclib_rates::rate_Cl35_to_p_S34_reaclib(rate_eval, tfactors);
        }

        {
            // Cl35_to_He4_P31_reaclib
            reaclib_rates::rate_Cl35_to_He4_P31_reaclib(rate_eval, tfactors);
        }

        {
            // Cl37_to_He4_P33_reaclib
            reaclib_rates::rate_Cl37_to_He4_P33_reaclib(rate_eval, tfactors);
        }

        {
            // Ar36_to_p_Cl35_reaclib
            reaclib_rates::rate_Ar36_to_p_Cl35_reaclib(rate_eval, tfactors);
        }

        {
            // Ar36_to_He4_S32_reaclib
            reaclib_rates::rate_Ar36_to_He4_S32_reaclib(rate_eval, tfactors);
        }

        {
            // Ar37_to_n_Ar36_reaclib
            reaclib_rates::rate_Ar37_to_n_Ar36_reaclib(rate_eval, tfactors);
        }

        {
            // Ar37_to_He4_S33_reaclib
            reaclib_rates::rate_Ar37_to_He4_S33_reaclib(rate_eval, tfactors);
        }

        {
            // Ar38_to_n_Ar37_reaclib
            reaclib_rates::rate_Ar38_to_n_Ar37_reaclib(rate_eval, tfactors);
        }

        {
            // Ar38_to_p_Cl37_reaclib
            reaclib_rates::rate_Ar38_to_p_Cl37_reaclib(rate_eval, tfactors);
        }

        {
            // Ar38_to_He4_S34_reaclib
            reaclib_rates::rate_Ar38_to_He4_S34_reaclib(rate_eval, tfactors);
        }

        {
            // K37_to_p_Ar36_reaclib
            reaclib_rates::rate_K37_to_p_Ar36_reaclib(rate_eval, tfactors);
        }

        {
            // K37_to_He4_Cl33_reaclib
            reaclib_rates::rate_K37_to_He4_Cl33_reaclib(rate_eval, tfactors);
        }

        {
            // K39_to_p_Ar38_reaclib
            reaclib_rates::rate_K39_to_p_Ar38_reaclib(rate_eval, tfactors);
        }

        {
            // K39_to_He4_Cl35_reaclib
            reaclib_rates::rate_K39_to_He4_Cl35_reaclib(rate_eval, tfactors);
        }

        {
            // Ca40_to_p_K39_reaclib
            reaclib_rates::rate_Ca40_to_p_K39_reaclib(rate_eval, tfactors);
        }

        {
            // Ca40_to_He4_Ar36_reaclib
            reaclib_rates::rate_Ca40_to_He4_Ar36_reaclib(rate_eval, tfactors);
        }

        {
            // Ca42_to_He4_Ar38_reaclib
            reaclib_rates::rate_Ca42_to_He4_Ar38_reaclib(rate_eval, tfactors);
        }

        {
            // Sc41_to_p_Ca40_reaclib
            reaclib_rates::rate_Sc41_to_p_Ca40_reaclib(rate_eval, tfactors);
        }

        {
            // Sc41_to_He4_K37_reaclib
            reaclib_rates::rate_Sc41_to_He4_K37_reaclib(rate_eval, tfactors);
        }

        {
            // Sc43_to_p_Ca42_reaclib
            reaclib_rates::rate_Sc43_to_p_Ca42_reaclib(rate_eval, tfactors);
        }

        {
            // Sc43_to_He4_K39_reaclib
            reaclib_rates::rate_Sc43_to_He4_K39_reaclib(rate_eval, tfactors);
        }

        {
            // Ti44_to_p_Sc43_reaclib
            reaclib_rates::rate_Ti44_to_p_Sc43_reaclib(rate_eval, tfactors);
        }

        {
            // Ti44_to_He4_Ca40_reaclib
            reaclib_rates::rate_Ti44_to_He4_Ca40_reaclib(rate_eval, tfactors);
        }

        {
            // Ti46_to_p_Sc45_reaclib
            reaclib_rates::rate_Ti46_to_p_Sc45_reaclib(rate_eval, tfactors);
        }

        {
            // Ti46_to_He4_Ca42_reaclib
            reaclib_rates::rate_Ti46_to_He4_Ca42_reaclib(rate_eval, tfactors);
        }

        {
            // V47_to_p_Ti46_reaclib
            reaclib_rates::rate_V47_to_p_Ti46_reaclib(rate_eval, tfactors);
        }

        {
            // V47_to_He4_Sc43_reaclib
            reaclib_rates::rate_V47_to_He4_Sc43_reaclib(rate_eval, tfactors);
        }

        {
            // V49_to_p_Ti48_reaclib
            reaclib_rates::rate_V49_to_p_Ti48_reaclib(rate_eval, tfactors);
        }

        {
            // V49_to_He4_Sc45_reaclib
            reaclib_rates::rate_V49_to_He4_Sc45_reaclib(rate_eval, tfactors);
        }

        {
            // Cr48_to_p_V47_reaclib
            reaclib_rates::rate_Cr48_to_p_V47_reaclib(rate_eval, tfactors);
        }

        {
            // Cr48_to_He4_Ti44_reaclib
            reaclib_rates::rate_Cr48_to_He4_Ti44_reaclib(rate_eval, tfactors);
        }

        {
            // Cr50_to_p_V49_reaclib
            reaclib_rates::rate_Cr50_to_p_V49_reaclib(rate_eval, tfactors);
        }

        {
            // Cr50_to_He4_Ti46_reaclib
            reaclib_rates::rate_Cr50_to_He4_Ti46_reaclib(rate_eval, tfactors);
        }

        {
            // Cr52_to_p_V51_reaclib
            reaclib_rates::rate_Cr52_to_p_V51_reaclib(rate_eval, tfactors);
        }

        {
            // Cr52_to_He4_Ti48_reaclib
            reaclib_rates::rate_Cr52_to_He4_Ti48_reaclib(rate_eval, tfactors);
        }

        {
            // Mn51_to_p_Cr50_reaclib
            reaclib_rates::rate_Mn51_to_p_Cr50_reaclib(rate_eval, tfactors);
        }

        {
            // Mn51_to_He4_V47_reaclib
            reaclib_rates::rate_Mn51_to_He4_V47_reaclib(rate_eval, tfactors);
        }

        {
            // Mn53_to_p_Cr52_reaclib
            reaclib_rates::rate_Mn53_to_p_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // Mn53_to_He4_V49_reaclib
            reaclib_rates::rate_Mn53_to_He4_V49_reaclib(rate_eval, tfactors);
        }

        {
            // Mn55_to_p_Cr54_reaclib
            reaclib_rates::rate_Mn55_to_p_Cr54_reaclib(rate_eval, tfactors);
        }

        {
            // Mn55_to_He4_V51_reaclib
            reaclib_rates::rate_Mn55_to_He4_V51_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_to_p_Mn51_reaclib
            reaclib_rates::rate_Fe52_to_p_Mn51_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_to_He4_Cr48_reaclib
            reaclib_rates::rate_Fe52_to_He4_Cr48_reaclib(rate_eval, tfactors);
        }

        {
            // Fe54_to_p_Mn53_reaclib
            reaclib_rates::rate_Fe54_to_p_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe54_to_He4_Cr50_reaclib
            reaclib_rates::rate_Fe54_to_He4_Cr50_reaclib(rate_eval, tfactors);
        }

        {
            // Fe56_to_p_Mn55_reaclib
            reaclib_rates::rate_Fe56_to_p_Mn55_reaclib(rate_eval, tfactors);
        }

        {
            // Fe56_to_He4_Cr52_reaclib
            reaclib_rates::rate_Fe56_to_He4_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // Fe58_to_He4_Cr54_reaclib
            reaclib_rates::rate_Fe58_to_He4_Cr54_reaclib(rate_eval, tfactors);
        }

        {
            // Co55_to_p_Fe54_reaclib
            reaclib_rates::rate_Co55_to_p_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Co55_to_He4_Mn51_reaclib
            reaclib_rates::rate_Co55_to_He4_Mn51_reaclib(rate_eval, tfactors);
        }

        {
            // Co56_to_n_Co55_reaclib
            reaclib_rates::rate_Co56_to_n_Co55_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_to_n_Co56_reaclib
            reaclib_rates::rate_Co57_to_n_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_to_p_Fe56_reaclib
            reaclib_rates::rate_Co57_to_p_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_to_He4_Mn53_reaclib
            reaclib_rates::rate_Co57_to_He4_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Co59_to_p_Fe58_reaclib
            reaclib_rates::rate_Co59_to_p_Fe58_reaclib(rate_eval, tfactors);
        }

        {
            // Co59_to_He4_Mn55_reaclib
            reaclib_rates::rate_Co59_to_He4_Mn55_reaclib(rate_eval, tfactors);
        }

        {
            // Ni56_to_p_Co55_reaclib
            reaclib_rates::rate_Ni56_to_p_Co55_reaclib(rate_eval, tfactors);
        }

        {
            // Ni56_to_He4_Fe52_reaclib
            reaclib_rates::rate_Ni56_to_He4_Fe52_reaclib(rate_eval, tfactors);
        }

        {
            // Ni58_to_p_Co57_reaclib
            reaclib_rates::rate_Ni58_to_p_Co57_reaclib(rate_eval, tfactors);
        }

        {
            // Ni58_to_He4_Fe54_reaclib
            reaclib_rates::rate_Ni58_to_He4_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Ni60_to_p_Co59_reaclib
            reaclib_rates::rate_Ni60_to_p_Co59_reaclib(rate_eval, tfactors);
        }

        {
            // Ni60_to_He4_Fe56_reaclib
            reaclib_rates::rate_Ni60_to_He4_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Cu59_to_p_Ni58_reaclib
            reaclib_rates::rate_Cu59_to_p_Ni58_reaclib(rate_eval, tfactors);
        }

        {
            // Cu59_to_He4_Co55_reaclib
            reaclib_rates::rate_Cu59_to_He4_Co55_reaclib(rate_eval, tfactors);
        }

        {
            // Zn60_to_p_Cu59_reaclib
            reaclib_rates::rate_Zn60_to_p_Cu59_reaclib(rate_eval, tfactors);
        }

        {
            // Zn60_to_He4_Ni56_reaclib
            reaclib_rates::rate_Zn60_to_He4_Ni56_reaclib(rate_eval, tfactors);
        }

        {
            // C12_to_3He4_reaclib
            reaclib_rates::rate_C12_to_3He4_reaclib(rate_eval, tfactors);
        }

        {
            // C12_He4_to_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C12);
            }
#endif
            reaclib_rates::rate_C12_He4_to_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_He4_to_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O16);
            }
#endif
            reaclib_rates::rate_O16_He4_to_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg23_n_to_Mg24_reaclib
            reaclib_rates::rate_Mg23_n_to_Mg24_reaclib(rate_eval, tfactors);
        }

        {
            // Mg23_He4_to_Si27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg23);
            }
#endif
            reaclib_rates::rate_Mg23_He4_to_Si27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_He4_to_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_He4_to_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_p_to_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            reaclib_rates::rate_Al27_p_to_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_He4_to_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al27);
            }
#endif
            reaclib_rates::rate_Al27_He4_to_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si27_n_to_Si28_reaclib
            reaclib_rates::rate_Si27_n_to_Si28_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_n_to_Si29_reaclib
            reaclib_rates::rate_Si28_n_to_Si29_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_p_to_P29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si28);
            }
#endif
            reaclib_rates::rate_Si28_p_to_P29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si28_He4_to_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            reaclib_rates::rate_Si28_He4_to_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si29_n_to_Si30_reaclib
            reaclib_rates::rate_Si29_n_to_Si30_reaclib(rate_eval, tfactors);
        }

        {
            // Si29_He4_to_S33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si29);
            }
#endif
            reaclib_rates::rate_Si29_He4_to_S33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si30_n_to_Si31_reaclib
            reaclib_rates::rate_Si30_n_to_Si31_reaclib(rate_eval, tfactors);
        }

        {
            // Si30_p_to_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si30);
            }
#endif
            reaclib_rates::rate_Si30_p_to_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si30_He4_to_S34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si30);
            }
#endif
            reaclib_rates::rate_Si30_He4_to_S34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si31_He4_to_S35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si31);
            }
#endif
            reaclib_rates::rate_Si31_He4_to_S35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P29_He4_to_Cl33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P29);
            }
#endif
            reaclib_rates::rate_P29_He4_to_Cl33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_p_to_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            reaclib_rates::rate_P31_p_to_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_He4_to_Cl35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P31);
            }
#endif
            reaclib_rates::rate_P31_He4_to_Cl35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P33_p_to_S34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P33);
            }
#endif
            reaclib_rates::rate_P33_p_to_S34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P33_He4_to_Cl37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P33);
            }
#endif
            reaclib_rates::rate_P33_He4_to_Cl37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S32_n_to_S33_reaclib
            reaclib_rates::rate_S32_n_to_S33_reaclib(rate_eval, tfactors);
        }

        {
            // S32_p_to_Cl33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S32);
            }
#endif
            reaclib_rates::rate_S32_p_to_Cl33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S32_He4_to_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
            }
#endif
            reaclib_rates::rate_S32_He4_to_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S33_n_to_S34_reaclib
            reaclib_rates::rate_S33_n_to_S34_reaclib(rate_eval, tfactors);
        }

        {
            // S33_He4_to_Ar37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S33);
            }
#endif
            reaclib_rates::rate_S33_He4_to_Ar37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S34_n_to_S35_reaclib
            reaclib_rates::rate_S34_n_to_S35_reaclib(rate_eval, tfactors);
        }

        {
            // S34_p_to_Cl35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S34);
            }
#endif
            reaclib_rates::rate_S34_p_to_Cl35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S34_He4_to_Ar38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S34);
            }
#endif
            reaclib_rates::rate_S34_He4_to_Ar38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl33_He4_to_K37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl33);
            }
#endif
            reaclib_rates::rate_Cl33_He4_to_K37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_p_to_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl35);
            }
#endif
            reaclib_rates::rate_Cl35_p_to_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_He4_to_K39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl35);
            }
#endif
            reaclib_rates::rate_Cl35_He4_to_K39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl37_p_to_Ar38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl37);
            }
#endif
            reaclib_rates::rate_Cl37_p_to_Ar38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar36_n_to_Ar37_reaclib
            reaclib_rates::rate_Ar36_n_to_Ar37_reaclib(rate_eval, tfactors);
        }

        {
            // Ar36_p_to_K37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar36);
            }
#endif
            reaclib_rates::rate_Ar36_p_to_K37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar36_He4_to_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
            }
#endif
            reaclib_rates::rate_Ar36_He4_to_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar37_n_to_Ar38_reaclib
            reaclib_rates::rate_Ar37_n_to_Ar38_reaclib(rate_eval, tfactors);
        }

        {
            // Ar38_p_to_K39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar38);
            }
#endif
            reaclib_rates::rate_Ar38_p_to_K39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar38_He4_to_Ca42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar38);
            }
#endif
            reaclib_rates::rate_Ar38_He4_to_Ca42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K37_He4_to_Sc41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K37);
            }
#endif
            reaclib_rates::rate_K37_He4_to_Sc41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K39_p_to_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K39);
            }
#endif
            reaclib_rates::rate_K39_p_to_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K39_He4_to_Sc43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K39);
            }
#endif
            reaclib_rates::rate_K39_He4_to_Sc43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_p_to_Sc41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca40);
            }
#endif
            reaclib_rates::rate_Ca40_p_to_Sc41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_He4_to_Ti44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
            }
#endif
            reaclib_rates::rate_Ca40_He4_to_Ti44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca42_p_to_Sc43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca42);
            }
#endif
            reaclib_rates::rate_Ca42_p_to_Sc43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca42_He4_to_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca42);
            }
#endif
            reaclib_rates::rate_Ca42_He4_to_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_p_to_Ti44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc43);
            }
#endif
            reaclib_rates::rate_Sc43_p_to_Ti44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_He4_to_V47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc43);
            }
#endif
            reaclib_rates::rate_Sc43_He4_to_V47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc45_p_to_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc45);
            }
#endif
            reaclib_rates::rate_Sc45_p_to_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc45_He4_to_V49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc45);
            }
#endif
            reaclib_rates::rate_Sc45_He4_to_V49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_He4_to_Cr48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti44);
            }
#endif
            reaclib_rates::rate_Ti44_He4_to_Cr48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti46_p_to_V47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti46);
            }
#endif
            reaclib_rates::rate_Ti46_p_to_V47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti46_He4_to_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti46);
            }
#endif
            reaclib_rates::rate_Ti46_He4_to_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti48_p_to_V49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti48);
            }
#endif
            reaclib_rates::rate_Ti48_p_to_V49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti48_He4_to_Cr52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti48);
            }
#endif
            reaclib_rates::rate_Ti48_He4_to_Cr52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V47_p_to_Cr48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V47);
            }
#endif
            reaclib_rates::rate_V47_p_to_Cr48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V47_He4_to_Mn51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V47);
            }
#endif
            reaclib_rates::rate_V47_He4_to_Mn51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_p_to_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V49);
            }
#endif
            reaclib_rates::rate_V49_p_to_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_He4_to_Mn53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V49);
            }
#endif
            reaclib_rates::rate_V49_He4_to_Mn53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V51_p_to_Cr52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V51);
            }
#endif
            reaclib_rates::rate_V51_p_to_Cr52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V51_He4_to_Mn55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V51);
            }
#endif
            reaclib_rates::rate_V51_He4_to_Mn55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr48_He4_to_Fe52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr48);
            }
#endif
            reaclib_rates::rate_Cr48_He4_to_Fe52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr50_p_to_Mn51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr50);
            }
#endif
            reaclib_rates::rate_Cr50_p_to_Mn51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr50_He4_to_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr50);
            }
#endif
            reaclib_rates::rate_Cr50_He4_to_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr52_p_to_Mn53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr52);
            }
#endif
            reaclib_rates::rate_Cr52_p_to_Mn53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr52_He4_to_Fe56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr52);
            }
#endif
            reaclib_rates::rate_Cr52_He4_to_Fe56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr54_p_to_Mn55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr54);
            }
#endif
            reaclib_rates::rate_Cr54_p_to_Mn55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr54_He4_to_Fe58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr54);
            }
#endif
            reaclib_rates::rate_Cr54_He4_to_Fe58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_p_to_Fe52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
            }
#endif
            reaclib_rates::rate_Mn51_p_to_Fe52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_He4_to_Co55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn51);
            }
#endif
            reaclib_rates::rate_Mn51_He4_to_Co55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_p_to_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn53);
            }
#endif
            reaclib_rates::rate_Mn53_p_to_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_He4_to_Co57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn53);
            }
#endif
            reaclib_rates::rate_Mn53_He4_to_Co57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn55_p_to_Fe56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn55);
            }
#endif
            reaclib_rates::rate_Mn55_p_to_Fe56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn55_He4_to_Co59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn55);
            }
#endif
            reaclib_rates::rate_Mn55_He4_to_Co59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe52_He4_to_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe52);
            }
#endif
            reaclib_rates::rate_Fe52_He4_to_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe54_p_to_Co55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
            }
#endif
            reaclib_rates::rate_Fe54_p_to_Co55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe54_He4_to_Ni58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
            }
#endif
            reaclib_rates::rate_Fe54_He4_to_Ni58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe56_p_to_Co57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            reaclib_rates::rate_Fe56_p_to_Co57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe56_He4_to_Ni60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe56);
            }
#endif
            reaclib_rates::rate_Fe56_He4_to_Ni60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe58_p_to_Co59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe58);
            }
#endif
            reaclib_rates::rate_Fe58_p_to_Co59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_n_to_Co56_reaclib
            reaclib_rates::rate_Co55_n_to_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Co55_p_to_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
            }
#endif
            reaclib_rates::rate_Co55_p_to_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_He4_to_Cu59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co55);
            }
#endif
            reaclib_rates::rate_Co55_He4_to_Cu59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co56_n_to_Co57_reaclib
            reaclib_rates::rate_Co56_n_to_Co57_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_p_to_Ni58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
            }
#endif
            reaclib_rates::rate_Co57_p_to_Ni58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co59_p_to_Ni60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co59);
            }
#endif
            reaclib_rates::rate_Co59_p_to_Ni60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni56_He4_to_Zn60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
            }
#endif
            reaclib_rates::rate_Ni56_He4_to_Zn60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni58_p_to_Cu59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
            }
#endif
            reaclib_rates::rate_Ni58_p_to_Cu59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu59_p_to_Zn60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu59);
            }
#endif
            reaclib_rates::rate_Cu59_p_to_Zn60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C12_C12_to_n_Mg23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            reaclib_rates::rate_C12_C12_to_n_Mg23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_C12_to_n_Si27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            reaclib_rates::rate_O16_C12_to_n_Si27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_C12_to_p_Al27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            reaclib_rates::rate_O16_C12_to_p_Al27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_C12_to_He4_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            reaclib_rates::rate_O16_C12_to_He4_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_O16_to_p_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_O16_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_O16_O16);
            }
#endif
            reaclib_rates::rate_O16_O16_to_p_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_O16_to_He4_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_O16_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_O16_O16);
            }
#endif
            reaclib_rates::rate_O16_O16_to_He4_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_He4_to_n_Mg23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_He4_to_n_Mg23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_C12_to_p_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_C12_to_p_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_C12_to_He4_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_C12_to_He4_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg23_n_to_He4_Ne20_reaclib
            reaclib_rates::rate_Mg23_n_to_He4_Ne20_reaclib(rate_eval, tfactors);
        }

        {
            // Mg23_n_to_C12_C12_reaclib
            reaclib_rates::rate_Mg23_n_to_C12_C12_reaclib(rate_eval, tfactors);
        }

        {
            // Mg24_He4_to_n_Si27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_He4_to_n_Si27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_He4_to_p_Al27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_He4_to_p_Al27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_He4_to_C12_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_He4_to_C12_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_p_to_n_Si27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            reaclib_rates::rate_Al27_p_to_n_Si27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_p_to_He4_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            reaclib_rates::rate_Al27_p_to_He4_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_p_to_C12_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            reaclib_rates::rate_Al27_p_to_C12_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_He4_to_p_Si30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al27);
            }
#endif
            reaclib_rates::rate_Al27_He4_to_p_Si30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si27_n_to_p_Al27_reaclib
            reaclib_rates::rate_Si27_n_to_p_Al27_reaclib(rate_eval, tfactors);
        }

        {
            // Si27_n_to_He4_Mg24_reaclib
            reaclib_rates::rate_Si27_n_to_He4_Mg24_reaclib(rate_eval, tfactors);
        }

        {
            // Si27_n_to_C12_O16_reaclib
            reaclib_rates::rate_Si27_n_to_C12_O16_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_He4_to_p_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            reaclib_rates::rate_Si28_He4_to_p_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si28_He4_to_C12_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            reaclib_rates::rate_Si28_He4_to_C12_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si28_He4_to_O16_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            reaclib_rates::rate_Si28_He4_to_O16_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si29_p_to_n_P29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si29);
            }
#endif
            reaclib_rates::rate_Si29_p_to_n_P29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si29_He4_to_n_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si29);
            }
#endif
            reaclib_rates::rate_Si29_He4_to_n_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si30_p_to_He4_Al27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si30);
            }
#endif
            reaclib_rates::rate_Si30_p_to_He4_Al27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si30_He4_to_n_S33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si30);
            }
#endif
            reaclib_rates::rate_Si30_He4_to_n_S33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si30_He4_to_p_P33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si30);
            }
#endif
            reaclib_rates::rate_Si30_He4_to_p_P33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si31_p_to_n_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si31);
            }
#endif
            reaclib_rates::rate_Si31_p_to_n_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si31_He4_to_n_S34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si31);
            }
#endif
            reaclib_rates::rate_Si31_He4_to_n_S34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P29_n_to_p_Si29_reaclib
            reaclib_rates::rate_P29_n_to_p_Si29_reaclib(rate_eval, tfactors);
        }

        {
            // P29_He4_to_p_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P29);
            }
#endif
            reaclib_rates::rate_P29_He4_to_p_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_n_to_p_Si31_reaclib
            reaclib_rates::rate_P31_n_to_p_Si31_reaclib(rate_eval, tfactors);
        }

        {
            // P31_p_to_He4_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            reaclib_rates::rate_P31_p_to_He4_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_p_to_C12_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            reaclib_rates::rate_P31_p_to_C12_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_p_to_O16_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            reaclib_rates::rate_P31_p_to_O16_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_He4_to_p_S34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P31);
            }
#endif
            reaclib_rates::rate_P31_He4_to_p_S34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P33_p_to_n_S33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P33);
            }
#endif
            reaclib_rates::rate_P33_p_to_n_S33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P33_p_to_He4_Si30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P33);
            }
#endif
            reaclib_rates::rate_P33_p_to_He4_Si30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S32_n_to_He4_Si29_reaclib
            reaclib_rates::rate_S32_n_to_He4_Si29_reaclib(rate_eval, tfactors);
        }

        {
            // S32_p_to_He4_P29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S32);
            }
#endif
            reaclib_rates::rate_S32_p_to_He4_P29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S32_He4_to_p_Cl35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
            }
#endif
            reaclib_rates::rate_S32_He4_to_p_Cl35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S33_n_to_p_P33_reaclib
            reaclib_rates::rate_S33_n_to_p_P33_reaclib(rate_eval, tfactors);
        }

        {
            // S33_n_to_He4_Si30_reaclib
            reaclib_rates::rate_S33_n_to_He4_Si30_reaclib(rate_eval, tfactors);
        }

        {
            // S33_p_to_n_Cl33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S33);
            }
#endif
            reaclib_rates::rate_S33_p_to_n_Cl33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S33_He4_to_n_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S33);
            }
#endif
            reaclib_rates::rate_S33_He4_to_n_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S34_n_to_He4_Si31_reaclib
            reaclib_rates::rate_S34_n_to_He4_Si31_reaclib(rate_eval, tfactors);
        }

        {
            // S34_p_to_He4_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S34);
            }
#endif
            reaclib_rates::rate_S34_p_to_He4_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S34_He4_to_n_Ar37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S34);
            }
#endif
            reaclib_rates::rate_S34_He4_to_n_Ar37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S34_He4_to_p_Cl37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S34);
            }
#endif
            reaclib_rates::rate_S34_He4_to_p_Cl37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S35_p_to_n_Cl35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S35);
            }
#endif
            reaclib_rates::rate_S35_p_to_n_Cl35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S35_He4_to_n_Ar38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S35);
            }
#endif
            reaclib_rates::rate_S35_He4_to_n_Ar38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl33_n_to_p_S33_reaclib
            reaclib_rates::rate_Cl33_n_to_p_S33_reaclib(rate_eval, tfactors);
        }

        {
            // Cl33_He4_to_p_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl33);
            }
#endif
            reaclib_rates::rate_Cl33_He4_to_p_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_n_to_p_S35_reaclib
            reaclib_rates::rate_Cl35_n_to_p_S35_reaclib(rate_eval, tfactors);
        }

        {
            // Cl35_p_to_He4_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl35);
            }
#endif
            reaclib_rates::rate_Cl35_p_to_He4_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_He4_to_p_Ar38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl35);
            }
#endif
            reaclib_rates::rate_Cl35_He4_to_p_Ar38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl37_p_to_n_Ar37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl37);
            }
#endif
            reaclib_rates::rate_Cl37_p_to_n_Ar37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl37_p_to_He4_S34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl37);
            }
#endif
            reaclib_rates::rate_Cl37_p_to_He4_S34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar36_n_to_He4_S33_reaclib
            reaclib_rates::rate_Ar36_n_to_He4_S33_reaclib(rate_eval, tfactors);
        }

        {
            // Ar36_p_to_He4_Cl33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar36);
            }
#endif
            reaclib_rates::rate_Ar36_p_to_He4_Cl33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar36_He4_to_p_K39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
            }
#endif
            reaclib_rates::rate_Ar36_He4_to_p_K39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar37_n_to_p_Cl37_reaclib
            reaclib_rates::rate_Ar37_n_to_p_Cl37_reaclib(rate_eval, tfactors);
        }

        {
            // Ar37_n_to_He4_S34_reaclib
            reaclib_rates::rate_Ar37_n_to_He4_S34_reaclib(rate_eval, tfactors);
        }

        {
            // Ar37_p_to_n_K37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar37);
            }
#endif
            reaclib_rates::rate_Ar37_p_to_n_K37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar37_He4_to_n_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar37);
            }
#endif
            reaclib_rates::rate_Ar37_He4_to_n_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar38_n_to_He4_S35_reaclib
            reaclib_rates::rate_Ar38_n_to_He4_S35_reaclib(rate_eval, tfactors);
        }

        {
            // Ar38_p_to_He4_Cl35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar38);
            }
#endif
            reaclib_rates::rate_Ar38_p_to_He4_Cl35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K37_n_to_p_Ar37_reaclib
            reaclib_rates::rate_K37_n_to_p_Ar37_reaclib(rate_eval, tfactors);
        }

        {
            // K37_He4_to_p_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K37);
            }
#endif
            reaclib_rates::rate_K37_He4_to_p_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K39_p_to_He4_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K39);
            }
#endif
            reaclib_rates::rate_K39_p_to_He4_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K39_He4_to_p_Ca42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K39);
            }
#endif
            reaclib_rates::rate_K39_He4_to_p_Ca42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_n_to_He4_Ar37_reaclib
            reaclib_rates::rate_Ca40_n_to_He4_Ar37_reaclib(rate_eval, tfactors);
        }

        {
            // Ca40_p_to_He4_K37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca40);
            }
#endif
            reaclib_rates::rate_Ca40_p_to_He4_K37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_He4_to_p_Sc43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
            }
#endif
            reaclib_rates::rate_Ca40_He4_to_p_Sc43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca42_p_to_He4_K39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca42);
            }
#endif
            reaclib_rates::rate_Ca42_p_to_He4_K39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca42_He4_to_p_Sc45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca42);
            }
#endif
            reaclib_rates::rate_Ca42_He4_to_p_Sc45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc41_He4_to_p_Ti44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc41);
            }
#endif
            reaclib_rates::rate_Sc41_He4_to_p_Ti44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_p_to_He4_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc43);
            }
#endif
            reaclib_rates::rate_Sc43_p_to_He4_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_He4_to_p_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc43);
            }
#endif
            reaclib_rates::rate_Sc43_He4_to_p_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc45_p_to_He4_Ca42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc45);
            }
#endif
            reaclib_rates::rate_Sc45_p_to_He4_Ca42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc45_He4_to_p_Ti48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc45);
            }
#endif
            reaclib_rates::rate_Sc45_He4_to_p_Ti48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_p_to_He4_Sc41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti44);
            }
#endif
            reaclib_rates::rate_Ti44_p_to_He4_Sc41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_He4_to_p_V47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti44);
            }
#endif
            reaclib_rates::rate_Ti44_He4_to_p_V47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti46_p_to_He4_Sc43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti46);
            }
#endif
            reaclib_rates::rate_Ti46_p_to_He4_Sc43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti46_He4_to_p_V49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti46);
            }
#endif
            reaclib_rates::rate_Ti46_He4_to_p_V49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti48_p_to_He4_Sc45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti48);
            }
#endif
            reaclib_rates::rate_Ti48_p_to_He4_Sc45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti48_He4_to_p_V51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti48);
            }
#endif
            reaclib_rates::rate_Ti48_He4_to_p_V51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V47_p_to_He4_Ti44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V47);
            }
#endif
            reaclib_rates::rate_V47_p_to_He4_Ti44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V47_He4_to_p_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V47);
            }
#endif
            reaclib_rates::rate_V47_He4_to_p_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_p_to_He4_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V49);
            }
#endif
            reaclib_rates::rate_V49_p_to_He4_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_He4_to_p_Cr52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V49);
            }
#endif
            reaclib_rates::rate_V49_He4_to_p_Cr52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V51_p_to_He4_Ti48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V51);
            }
#endif
            reaclib_rates::rate_V51_p_to_He4_Ti48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V51_He4_to_p_Cr54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V51);
            }
#endif
            reaclib_rates::rate_V51_He4_to_p_Cr54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr48_He4_to_p_Mn51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr48);
            }
#endif
            reaclib_rates::rate_Cr48_He4_to_p_Mn51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr50_p_to_He4_V47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr50);
            }
#endif
            reaclib_rates::rate_Cr50_p_to_He4_V47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr50_He4_to_p_Mn53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr50);
            }
#endif
            reaclib_rates::rate_Cr50_He4_to_p_Mn53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr52_p_to_He4_V49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr52);
            }
#endif
            reaclib_rates::rate_Cr52_p_to_He4_V49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr52_He4_to_p_Mn55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr52);
            }
#endif
            reaclib_rates::rate_Cr52_He4_to_p_Mn55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr54_p_to_He4_V51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr54);
            }
#endif
            reaclib_rates::rate_Cr54_p_to_He4_V51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_p_to_He4_Cr48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
            }
#endif
            reaclib_rates::rate_Mn51_p_to_He4_Cr48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_He4_to_p_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn51);
            }
#endif
            reaclib_rates::rate_Mn51_He4_to_p_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_p_to_He4_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn53);
            }
#endif
            reaclib_rates::rate_Mn53_p_to_He4_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_He4_to_n_Co56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn53);
            }
#endif
            reaclib_rates::rate_Mn53_He4_to_n_Co56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_He4_to_p_Fe56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn53);
            }
#endif
            reaclib_rates::rate_Mn53_He4_to_p_Fe56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn55_p_to_He4_Cr52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn55);
            }
#endif
            reaclib_rates::rate_Mn55_p_to_He4_Cr52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn55_He4_to_p_Fe58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn55);
            }
#endif
            reaclib_rates::rate_Mn55_He4_to_p_Fe58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe52_He4_to_p_Co55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe52);
            }
#endif
            reaclib_rates::rate_Fe52_He4_to_p_Co55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe54_p_to_He4_Mn51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
            }
#endif
            reaclib_rates::rate_Fe54_p_to_He4_Mn51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe54_He4_to_p_Co57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
            }
#endif
            reaclib_rates::rate_Fe54_He4_to_p_Co57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe56_p_to_n_Co56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            reaclib_rates::rate_Fe56_p_to_n_Co56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe56_p_to_He4_Mn53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            reaclib_rates::rate_Fe56_p_to_He4_Mn53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe56_He4_to_p_Co59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe56);
            }
#endif
            reaclib_rates::rate_Fe56_He4_to_p_Co59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe58_p_to_He4_Mn55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe58);
            }
#endif
            reaclib_rates::rate_Fe58_p_to_He4_Mn55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_p_to_He4_Fe52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
            }
#endif
            reaclib_rates::rate_Co55_p_to_He4_Fe52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_He4_to_p_Ni58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co55);
            }
#endif
            reaclib_rates::rate_Co55_He4_to_p_Ni58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co56_n_to_p_Fe56_reaclib
            reaclib_rates::rate_Co56_n_to_p_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Co56_n_to_He4_Mn53_reaclib
            reaclib_rates::rate_Co56_n_to_He4_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Co56_p_to_n_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
            }
#endif
            reaclib_rates::rate_Co56_p_to_n_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co56_He4_to_n_Cu59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co56);
            }
#endif
            reaclib_rates::rate_Co56_He4_to_n_Cu59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co57_p_to_He4_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
            }
#endif
            reaclib_rates::rate_Co57_p_to_He4_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co57_He4_to_p_Ni60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co57);
            }
#endif
            reaclib_rates::rate_Co57_He4_to_p_Ni60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co59_p_to_He4_Fe56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co59);
            }
#endif
            reaclib_rates::rate_Co59_p_to_He4_Fe56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni56_n_to_p_Co56_reaclib
            reaclib_rates::rate_Ni56_n_to_p_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Ni56_He4_to_p_Cu59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
            }
#endif
            reaclib_rates::rate_Ni56_He4_to_p_Cu59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni58_p_to_He4_Co55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
            }
#endif
            reaclib_rates::rate_Ni58_p_to_He4_Co55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni60_p_to_He4_Co57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni60);
            }
#endif
            reaclib_rates::rate_Ni60_p_to_He4_Co57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu59_n_to_He4_Co56_reaclib
            reaclib_rates::rate_Cu59_n_to_He4_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Cu59_p_to_He4_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu59);
            }
#endif
            reaclib_rates::rate_Cu59_p_to_He4_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He4_He4_He4_to_C12_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_He4_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_He4_Be8);
            }
#endif
            reaclib_rates::rate_He4_He4_He4_to_C12_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni56_n_to_Ni57_reaclib
            reaclib_rates::rate_Ni56_n_to_Ni57_reaclib(rate_eval, tfactors);
        }

        {
            // Ni57_to_n_Ni56_reaclib
            reaclib_rates::rate_Ni57_to_n_Ni56_reaclib(rate_eval, tfactors);
        }

        {
            // Ni57_n_to_Ni58_reaclib
            reaclib_rates::rate_Ni57_n_to_Ni58_reaclib(rate_eval, tfactors);
        }

        {
            // Ni58_to_n_Ni57_reaclib
            reaclib_rates::rate_Ni58_to_n_Ni57_reaclib(rate_eval, tfactors);
        }

        {
            // Ni58_n_to_Ni59_reaclib
            reaclib_rates::rate_Ni58_n_to_Ni59_reaclib(rate_eval, tfactors);
        }

        {
            // Ni59_to_n_Ni58_reaclib
            reaclib_rates::rate_Ni59_to_n_Ni58_reaclib(rate_eval, tfactors);
        }

        {
            // Ni59_n_to_Ni60_reaclib
            reaclib_rates::rate_Ni59_n_to_Ni60_reaclib(rate_eval, tfactors);
        }

        {
            // Ni60_to_n_Ni59_reaclib
            reaclib_rates::rate_Ni60_to_n_Ni59_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_n_to_Co58_reaclib
            reaclib_rates::rate_Co57_n_to_Co58_reaclib(rate_eval, tfactors);
        }

        {
            // Co58_to_n_Co57_reaclib
            reaclib_rates::rate_Co58_to_n_Co57_reaclib(rate_eval, tfactors);
        }

        {
            // Co58_n_to_Co59_reaclib
            reaclib_rates::rate_Co58_n_to_Co59_reaclib(rate_eval, tfactors);
        }

        {
            // Co59_to_n_Co58_reaclib
            reaclib_rates::rate_Co59_to_n_Co58_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_n_to_Fe53_reaclib
            reaclib_rates::rate_Fe52_n_to_Fe53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_to_n_Fe52_reaclib
            reaclib_rates::rate_Fe53_to_n_Fe52_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_n_to_Fe54_reaclib
            reaclib_rates::rate_Fe53_n_to_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Fe54_to_n_Fe53_reaclib
            reaclib_rates::rate_Fe54_to_n_Fe53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe54_n_to_Fe55_reaclib
            reaclib_rates::rate_Fe54_n_to_Fe55_reaclib(rate_eval, tfactors);
        }

        {
            // Fe55_to_n_Fe54_reaclib
            reaclib_rates::rate_Fe55_to_n_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Fe55_n_to_Fe56_reaclib
            reaclib_rates::rate_Fe55_n_to_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Fe56_to_n_Fe55_reaclib
            reaclib_rates::rate_Fe56_to_n_Fe55_reaclib(rate_eval, tfactors);
        }

        {
            // Fe56_n_to_Fe57_reaclib
            reaclib_rates::rate_Fe56_n_to_Fe57_reaclib(rate_eval, tfactors);
        }

        {
            // Fe57_to_n_Fe56_reaclib
            reaclib_rates::rate_Fe57_to_n_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Fe57_n_to_Fe58_reaclib
            reaclib_rates::rate_Fe57_n_to_Fe58_reaclib(rate_eval, tfactors);
        }

        {
            // Fe58_to_n_Fe57_reaclib
            reaclib_rates::rate_Fe58_to_n_Fe57_reaclib(rate_eval, tfactors);
        }

        {
            // Mn51_n_to_Mn52_reaclib
            reaclib_rates::rate_Mn51_n_to_Mn52_reaclib(rate_eval, tfactors);
        }

        {
            // Mn52_to_n_Mn51_reaclib
            reaclib_rates::rate_Mn52_to_n_Mn51_reaclib(rate_eval, tfactors);
        }

        {
            // Mn52_n_to_Mn53_reaclib
            reaclib_rates::rate_Mn52_n_to_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Mn53_to_n_Mn52_reaclib
            reaclib_rates::rate_Mn53_to_n_Mn52_reaclib(rate_eval, tfactors);
        }

        {
            // Mn53_n_to_Mn54_reaclib
            reaclib_rates::rate_Mn53_n_to_Mn54_reaclib(rate_eval, tfactors);
        }

        {
            // Mn54_to_n_Mn53_reaclib
            reaclib_rates::rate_Mn54_to_n_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Mn54_n_to_Mn55_reaclib
            reaclib_rates::rate_Mn54_n_to_Mn55_reaclib(rate_eval, tfactors);
        }

        {
            // Mn55_to_n_Mn54_reaclib
            reaclib_rates::rate_Mn55_to_n_Mn54_reaclib(rate_eval, tfactors);
        }

        {
            // Cr48_n_to_Cr49_reaclib
            reaclib_rates::rate_Cr48_n_to_Cr49_reaclib(rate_eval, tfactors);
        }

        {
            // Cr49_to_n_Cr48_reaclib
            reaclib_rates::rate_Cr49_to_n_Cr48_reaclib(rate_eval, tfactors);
        }

        {
            // Cr49_n_to_Cr50_reaclib
            reaclib_rates::rate_Cr49_n_to_Cr50_reaclib(rate_eval, tfactors);
        }

        {
            // Cr50_to_n_Cr49_reaclib
            reaclib_rates::rate_Cr50_to_n_Cr49_reaclib(rate_eval, tfactors);
        }

        {
            // Cr50_n_to_Cr51_reaclib
            reaclib_rates::rate_Cr50_n_to_Cr51_reaclib(rate_eval, tfactors);
        }

        {
            // Cr51_to_n_Cr50_reaclib
            reaclib_rates::rate_Cr51_to_n_Cr50_reaclib(rate_eval, tfactors);
        }

        {
            // Cr51_n_to_Cr52_reaclib
            reaclib_rates::rate_Cr51_n_to_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // Cr52_to_n_Cr51_reaclib
            reaclib_rates::rate_Cr52_to_n_Cr51_reaclib(rate_eval, tfactors);
        }

        {
            // Cr52_n_to_Cr53_reaclib
            reaclib_rates::rate_Cr52_n_to_Cr53_reaclib(rate_eval, tfactors);
        }

        {
            // Cr53_to_n_Cr52_reaclib
            reaclib_rates::rate_Cr53_to_n_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // Cr53_n_to_Cr54_reaclib
            reaclib_rates::rate_Cr53_n_to_Cr54_reaclib(rate_eval, tfactors);
        }

        {
            // Cr54_to_n_Cr53_reaclib
            reaclib_rates::rate_Cr54_to_n_Cr53_reaclib(rate_eval, tfactors);
        }

        {
            // V47_n_to_V48_reaclib
            reaclib_rates::rate_V47_n_to_V48_reaclib(rate_eval, tfactors);
        }

        {
            // V48_to_n_V47_reaclib
            reaclib_rates::rate_V48_to_n_V47_reaclib(rate_eval, tfactors);
        }

        {
            // V48_n_to_V49_reaclib
            reaclib_rates::rate_V48_n_to_V49_reaclib(rate_eval, tfactors);
        }

        {
            // V49_to_n_V48_reaclib
            reaclib_rates::rate_V49_to_n_V48_reaclib(rate_eval, tfactors);
        }

        {
            // V49_n_to_V50_reaclib
            reaclib_rates::rate_V49_n_to_V50_reaclib(rate_eval, tfactors);
        }

        {
            // V50_to_n_V49_reaclib
            reaclib_rates::rate_V50_to_n_V49_reaclib(rate_eval, tfactors);
        }

        {
            // V50_n_to_V51_reaclib
            reaclib_rates::rate_V50_n_to_V51_reaclib(rate_eval, tfactors);
        }

        {
            // V51_to_n_V50_reaclib
            reaclib_rates::rate_V51_to_n_V50_reaclib(rate_eval, tfactors);
        }

        {
            // Ti44_n_to_Ti45_reaclib
            reaclib_rates::rate_Ti44_n_to_Ti45_reaclib(rate_eval, tfactors);
        }

        {
            // Ti45_to_n_Ti44_reaclib
            reaclib_rates::rate_Ti45_to_n_Ti44_reaclib(rate_eval, tfactors);
        }

        {
            // Ti45_n_to_Ti46_reaclib
            reaclib_rates::rate_Ti45_n_to_Ti46_reaclib(rate_eval, tfactors);
        }

        {
            // Ti46_to_n_Ti45_reaclib
            reaclib_rates::rate_Ti46_to_n_Ti45_reaclib(rate_eval, tfactors);
        }

        {
            // Ti46_n_to_Ti47_reaclib
            reaclib_rates::rate_Ti46_n_to_Ti47_reaclib(rate_eval, tfactors);
        }

        {
            // Ti47_to_n_Ti46_reaclib
            reaclib_rates::rate_Ti47_to_n_Ti46_reaclib(rate_eval, tfactors);
        }

        {
            // Ti47_n_to_Ti48_reaclib
            reaclib_rates::rate_Ti47_n_to_Ti48_reaclib(rate_eval, tfactors);
        }

        {
            // Ti48_to_n_Ti47_reaclib
            reaclib_rates::rate_Ti48_to_n_Ti47_reaclib(rate_eval, tfactors);
        }

        {
            // Sc41_n_to_Sc42_reaclib
            reaclib_rates::rate_Sc41_n_to_Sc42_reaclib(rate_eval, tfactors);
        }

        {
            // Sc42_to_n_Sc41_reaclib
            reaclib_rates::rate_Sc42_to_n_Sc41_reaclib(rate_eval, tfactors);
        }

        {
            // Sc42_n_to_Sc43_reaclib
            reaclib_rates::rate_Sc42_n_to_Sc43_reaclib(rate_eval, tfactors);
        }

        {
            // Sc43_to_n_Sc42_reaclib
            reaclib_rates::rate_Sc43_to_n_Sc42_reaclib(rate_eval, tfactors);
        }

        {
            // Sc43_n_to_Sc44_reaclib
            reaclib_rates::rate_Sc43_n_to_Sc44_reaclib(rate_eval, tfactors);
        }

        {
            // Sc44_to_n_Sc43_reaclib
            reaclib_rates::rate_Sc44_to_n_Sc43_reaclib(rate_eval, tfactors);
        }

        {
            // Sc44_n_to_Sc45_reaclib
            reaclib_rates::rate_Sc44_n_to_Sc45_reaclib(rate_eval, tfactors);
        }

        {
            // Sc45_to_n_Sc44_reaclib
            reaclib_rates::rate_Sc45_to_n_Sc44_reaclib(rate_eval, tfactors);
        }

        {
            // Ca40_n_to_Ca41_reaclib
            reaclib_rates::rate_Ca40_n_to_Ca41_reaclib(rate_eval, tfactors);
        }

        {
            // Ca41_to_n_Ca40_reaclib
            reaclib_rates::rate_Ca41_to_n_Ca40_reaclib(rate_eval, tfactors);
        }

        {
            // Ca41_n_to_Ca42_reaclib
            reaclib_rates::rate_Ca41_n_to_Ca42_reaclib(rate_eval, tfactors);
        }

        {
            // Ca42_to_n_Ca41_reaclib
            reaclib_rates::rate_Ca42_to_n_Ca41_reaclib(rate_eval, tfactors);
        }

        {
            // C12_C12_to_p_Na23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            reaclib_rates::rate_C12_C12_to_p_Na23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na23_p_to_C12_C12_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            reaclib_rates::rate_Na23_p_to_C12_C12_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na23_p_to_He4_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            reaclib_rates::rate_Na23_p_to_He4_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na23_p_to_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            reaclib_rates::rate_Na23_p_to_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_to_p_Na23_reaclib
            reaclib_rates::rate_Mg24_to_p_Na23_reaclib(rate_eval, tfactors);
        }

        {
            // C12_C12_to_He4_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            reaclib_rates::rate_C12_C12_to_He4_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_He4_to_C12_C12_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_He4_to_C12_C12_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_He4_to_p_Na23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_He4_to_p_Na23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_He4_to_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_He4_to_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_to_He4_Ne20_reaclib
            reaclib_rates::rate_Mg24_to_He4_Ne20_reaclib(rate_eval, tfactors);
        }

        {
            // Al27_to_Si27_beta_neg_weaktab
            tabular_weak_rates::rate_Al27_to_Si27_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si27_to_Al27_electron_capture_weaktab
            tabular_weak_rates::rate_Si27_to_Al27_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P29_to_Si29_electron_capture_weaktab
            tabular_weak_rates::rate_P29_to_Si29_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si29_to_P29_beta_neg_weaktab
            tabular_weak_rates::rate_Si29_to_P29_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P31_to_Si31_electron_capture_weaktab
            tabular_weak_rates::rate_P31_to_Si31_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si31_to_P31_beta_neg_weaktab
            tabular_weak_rates::rate_Si31_to_P31_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl33_to_S33_electron_capture_weaktab
            tabular_weak_rates::rate_Cl33_to_S33_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P33_to_S33_beta_neg_weaktab
            tabular_weak_rates::rate_P33_to_S33_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S33_to_Cl33_beta_neg_weaktab
            tabular_weak_rates::rate_S33_to_Cl33_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S33_to_P33_electron_capture_weaktab
            tabular_weak_rates::rate_S33_to_P33_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl35_to_S35_electron_capture_weaktab
            tabular_weak_rates::rate_Cl35_to_S35_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S35_to_Cl35_beta_neg_weaktab
            tabular_weak_rates::rate_S35_to_Cl35_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar37_to_Cl37_electron_capture_weaktab
            tabular_weak_rates::rate_Ar37_to_Cl37_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar37_to_K37_beta_neg_weaktab
            tabular_weak_rates::rate_Ar37_to_K37_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl37_to_Ar37_beta_neg_weaktab
            tabular_weak_rates::rate_Cl37_to_Ar37_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K37_to_Ar37_electron_capture_weaktab
            tabular_weak_rates::rate_K37_to_Ar37_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co56_to_Fe56_electron_capture_weaktab
            tabular_weak_rates::rate_Co56_to_Fe56_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co56_to_Ni56_beta_neg_weaktab
            tabular_weak_rates::rate_Co56_to_Ni56_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe56_to_Co56_beta_neg_weaktab
            tabular_weak_rates::rate_Fe56_to_Co56_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni56_to_Co56_electron_capture_weaktab
            tabular_weak_rates::rate_Ni56_to_Co56_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // n_to_p_beta_neg_weaktab
            tabular_weak_rates::rate_n_to_p_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // p_to_n_electron_capture_weaktab
            tabular_weak_rates::rate_p_to_n_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni56_n_n_to_Ni58_approx
            approximate_rates::rate_Ni56_n_n_to_Ni58_approx(rate_eval, rho, Y);
        }

        {
            // Ni58_to_Ni56_n_n_approx
            approximate_rates::rate_Ni58_to_Ni56_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ni58_n_n_to_Ni60_approx
            approximate_rates::rate_Ni58_n_n_to_Ni60_approx(rate_eval, rho, Y);
        }

        {
            // Ni60_to_Ni58_n_n_approx
            approximate_rates::rate_Ni60_to_Ni58_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Co57_n_n_to_Co59_approx
            approximate_rates::rate_Co57_n_n_to_Co59_approx(rate_eval, rho, Y);
        }

        {
            // Co59_to_Co57_n_n_approx
            approximate_rates::rate_Co59_to_Co57_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Fe52_n_n_to_Fe54_approx
            approximate_rates::rate_Fe52_n_n_to_Fe54_approx(rate_eval, rho, Y);
        }

        {
            // Fe54_to_Fe52_n_n_approx
            approximate_rates::rate_Fe54_to_Fe52_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Fe54_n_n_to_Fe56_approx
            approximate_rates::rate_Fe54_n_n_to_Fe56_approx(rate_eval, rho, Y);
        }

        {
            // Fe56_to_Fe54_n_n_approx
            approximate_rates::rate_Fe56_to_Fe54_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Fe56_n_n_to_Fe58_approx
            approximate_rates::rate_Fe56_n_n_to_Fe58_approx(rate_eval, rho, Y);
        }

        {
            // Fe58_to_Fe56_n_n_approx
            approximate_rates::rate_Fe58_to_Fe56_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Mn51_n_n_to_Mn53_approx
            approximate_rates::rate_Mn51_n_n_to_Mn53_approx(rate_eval, rho, Y);
        }

        {
            // Mn53_to_Mn51_n_n_approx
            approximate_rates::rate_Mn53_to_Mn51_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Mn53_n_n_to_Mn55_approx
            approximate_rates::rate_Mn53_n_n_to_Mn55_approx(rate_eval, rho, Y);
        }

        {
            // Mn55_to_Mn53_n_n_approx
            approximate_rates::rate_Mn55_to_Mn53_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Cr48_n_n_to_Cr50_approx
            approximate_rates::rate_Cr48_n_n_to_Cr50_approx(rate_eval, rho, Y);
        }

        {
            // Cr50_to_Cr48_n_n_approx
            approximate_rates::rate_Cr50_to_Cr48_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Cr50_n_n_to_Cr52_approx
            approximate_rates::rate_Cr50_n_n_to_Cr52_approx(rate_eval, rho, Y);
        }

        {
            // Cr52_to_Cr50_n_n_approx
            approximate_rates::rate_Cr52_to_Cr50_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Cr52_n_n_to_Cr54_approx
            approximate_rates::rate_Cr52_n_n_to_Cr54_approx(rate_eval, rho, Y);
        }

        {
            // Cr54_to_Cr52_n_n_approx
            approximate_rates::rate_Cr54_to_Cr52_n_n_approx(rate_eval, rho, Y);
        }

        {
            // V47_n_n_to_V49_approx
            approximate_rates::rate_V47_n_n_to_V49_approx(rate_eval, rho, Y);
        }

        {
            // V49_to_V47_n_n_approx
            approximate_rates::rate_V49_to_V47_n_n_approx(rate_eval, rho, Y);
        }

        {
            // V49_n_n_to_V51_approx
            approximate_rates::rate_V49_n_n_to_V51_approx(rate_eval, rho, Y);
        }

        {
            // V51_to_V49_n_n_approx
            approximate_rates::rate_V51_to_V49_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ti44_n_n_to_Ti46_approx
            approximate_rates::rate_Ti44_n_n_to_Ti46_approx(rate_eval, rho, Y);
        }

        {
            // Ti46_to_Ti44_n_n_approx
            approximate_rates::rate_Ti46_to_Ti44_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ti46_n_n_to_Ti48_approx
            approximate_rates::rate_Ti46_n_n_to_Ti48_approx(rate_eval, rho, Y);
        }

        {
            // Ti48_to_Ti46_n_n_approx
            approximate_rates::rate_Ti48_to_Ti46_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Sc41_n_n_to_Sc43_approx
            approximate_rates::rate_Sc41_n_n_to_Sc43_approx(rate_eval, rho, Y);
        }

        {
            // Sc43_to_Sc41_n_n_approx
            approximate_rates::rate_Sc43_to_Sc41_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Sc43_n_n_to_Sc45_approx
            approximate_rates::rate_Sc43_n_n_to_Sc45_approx(rate_eval, rho, Y);
        }

        {
            // Sc45_to_Sc43_n_n_approx
            approximate_rates::rate_Sc45_to_Sc43_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ca40_n_n_to_Ca42_approx
            approximate_rates::rate_Ca40_n_n_to_Ca42_approx(rate_eval, rho, Y);
        }

        {
            // Ca42_to_Ca40_n_n_approx
            approximate_rates::rate_Ca42_to_Ca40_n_n_approx(rate_eval, rho, Y);
        }

        {
            // C12_C12_to_Mg24_approx
            approximate_rates::rate_C12_C12_to_Mg24_approx(rate_eval);
        }

        {
            // Mg24_to_C12_C12_approx
            approximate_rates::rate_Mg24_to_C12_C12_approx(rate_eval);
        }

        {
            // C12_C12_to_Ne20_He4_approx
            approximate_rates::rate_C12_C12_to_Ne20_He4_approx(rate_eval);
        }

        {
            // Ne20_He4_to_C12_C12_approx
            approximate_rates::rate_Ne20_He4_to_C12_C12_approx(rate_eval);
        }

        {
            // Ne20_He4_to_Mg24_approx
            approximate_rates::rate_Ne20_He4_to_Mg24_approx(rate_eval);
        }

        {
            // Mg24_to_Ne20_He4_approx
            approximate_rates::rate_Mg24_to_Ne20_He4_approx(rate_eval);
        }


    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_weak_rates([[maybe_unused]] const tf_t& tfactors,
                    [[maybe_unused]] const amrex::Real temp,
                    [[maybe_unused]] const amrex::Real rho,
                    [[maybe_unused]] const amrex::Real rhoy,
                    [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                    [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate{};
        [[maybe_unused]] amrex::Real drate_dT{};

        [[maybe_unused]] const amrex::Real log_temp = std::log10(temp);
        [[maybe_unused]] const amrex::Real log_rhoy = std::log10(rhoy);

        // we are not doing derivatives
        static_assert(std::is_same_v<T, rate_t>);

        rate_eval.enuc_weak = 0.0_rt;

        // rates are filled in a topologically-sorted order

        {
            // Al27_to_Si27_beta_neg_weaktab
            tabular_weak_rates::rate_Al27_to_Si27_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si27_to_Al27_electron_capture_weaktab
            tabular_weak_rates::rate_Si27_to_Al27_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P29_to_Si29_electron_capture_weaktab
            tabular_weak_rates::rate_P29_to_Si29_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si29_to_P29_beta_neg_weaktab
            tabular_weak_rates::rate_Si29_to_P29_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P31_to_Si31_electron_capture_weaktab
            tabular_weak_rates::rate_P31_to_Si31_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si31_to_P31_beta_neg_weaktab
            tabular_weak_rates::rate_Si31_to_P31_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl33_to_S33_electron_capture_weaktab
            tabular_weak_rates::rate_Cl33_to_S33_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P33_to_S33_beta_neg_weaktab
            tabular_weak_rates::rate_P33_to_S33_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S33_to_Cl33_beta_neg_weaktab
            tabular_weak_rates::rate_S33_to_Cl33_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S33_to_P33_electron_capture_weaktab
            tabular_weak_rates::rate_S33_to_P33_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl35_to_S35_electron_capture_weaktab
            tabular_weak_rates::rate_Cl35_to_S35_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S35_to_Cl35_beta_neg_weaktab
            tabular_weak_rates::rate_S35_to_Cl35_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar37_to_Cl37_electron_capture_weaktab
            tabular_weak_rates::rate_Ar37_to_Cl37_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar37_to_K37_beta_neg_weaktab
            tabular_weak_rates::rate_Ar37_to_K37_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl37_to_Ar37_beta_neg_weaktab
            tabular_weak_rates::rate_Cl37_to_Ar37_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K37_to_Ar37_electron_capture_weaktab
            tabular_weak_rates::rate_K37_to_Ar37_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co56_to_Fe56_electron_capture_weaktab
            tabular_weak_rates::rate_Co56_to_Fe56_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co56_to_Ni56_beta_neg_weaktab
            tabular_weak_rates::rate_Co56_to_Ni56_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe56_to_Co56_beta_neg_weaktab
            tabular_weak_rates::rate_Fe56_to_Co56_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni56_to_Co56_electron_capture_weaktab
            tabular_weak_rates::rate_Ni56_to_Co56_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // n_to_p_beta_neg_weaktab
            tabular_weak_rates::rate_n_to_p_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // p_to_n_electron_capture_weaktab
            tabular_weak_rates::rate_p_to_n_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }


    }

}

#endif
