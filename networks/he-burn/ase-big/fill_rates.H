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

        part_fun::pf_cache_t pf_cache{};

        pf_cache.index_temp_array_1 = interp_net::find_index(tfactors.T9, part_fun::temp_array_1);

        {
            // C14_to_N14_beta_neg_reaclib
            reaclib_rates::rate_C14_to_N14_beta_neg_reaclib(rate_eval, tfactors);
        }

        {
            // N13_to_C13_beta_pos_reaclib
            reaclib_rates::rate_N13_to_C13_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // C12_p_to_N13_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C12);
            }
#endif
            reaclib_rates::rate_C12_p_to_N13_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // C13_p_to_N14_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_C13);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C13);
            }
#endif
            reaclib_rates::rate_C13_p_to_N14_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C14_p_to_N15_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_C14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C14);
            }
#endif
            reaclib_rates::rate_C14_p_to_N15_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C14_He4_to_O18_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C14);
            }
#endif
            reaclib_rates::rate_C14_He4_to_O18_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N14_He4_to_F18_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N14);
            }
#endif
            reaclib_rates::rate_N14_He4_to_F18_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N15_p_to_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_N15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N15);
            }
#endif
            reaclib_rates::rate_N15_p_to_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N15_He4_to_F19_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N15);
            }
#endif
            reaclib_rates::rate_N15_He4_to_F19_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // O17_p_to_F18_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O17);
            }
#endif
            reaclib_rates::rate_O17_p_to_F18_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O17_He4_to_Ne21_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O17);
            }
#endif
            reaclib_rates::rate_O17_He4_to_Ne21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O18_p_to_F19_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O18);
            }
#endif
            reaclib_rates::rate_O18_p_to_F19_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O18_He4_to_Ne22_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O18);
            }
#endif
            reaclib_rates::rate_O18_He4_to_Ne22_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F18_He4_to_Na22_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F18);
            }
#endif
            reaclib_rates::rate_F18_He4_to_Na22_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F19_p_to_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F19);
            }
#endif
            reaclib_rates::rate_F19_p_to_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F19_He4_to_Na23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F19);
            }
#endif
            reaclib_rates::rate_F19_He4_to_Na23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ne21_p_to_Na22_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne21);
            }
#endif
            reaclib_rates::rate_Ne21_p_to_Na22_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne21_He4_to_Mg25_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne21);
            }
#endif
            reaclib_rates::rate_Ne21_He4_to_Mg25_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne22_p_to_Na23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne22);
            }
#endif
            reaclib_rates::rate_Ne22_p_to_Na23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na22_He4_to_Al26_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na22);
            }
#endif
            reaclib_rates::rate_Na22_He4_to_Al26_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Na23_He4_to_Al27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na23);
            }
#endif
            reaclib_rates::rate_Na23_He4_to_Al27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Mg25_p_to_Al26_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg25);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg25);
            }
#endif
            reaclib_rates::rate_Mg25_p_to_Al26_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg25_He4_to_Si29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg25);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg25);
            }
#endif
            reaclib_rates::rate_Mg25_He4_to_Si29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al26_He4_to_P30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al26);
            }
#endif
            reaclib_rates::rate_Al26_He4_to_P30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Si29_p_to_P30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si29);
            }
#endif
            reaclib_rates::rate_Si29_p_to_P30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // P30_He4_to_Cl34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P30);
            }
#endif
            reaclib_rates::rate_P30_He4_to_Cl34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // S33_p_to_Cl34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S33);
            }
#endif
            reaclib_rates::rate_S33_p_to_Cl34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Cl34_He4_to_K38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl34);
            }
#endif
            reaclib_rates::rate_Cl34_He4_to_K38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ar37_p_to_K38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar37);
            }
#endif
            reaclib_rates::rate_Ar37_p_to_K38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar37_He4_to_Ca41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar37);
            }
#endif
            reaclib_rates::rate_Ar37_He4_to_Ca41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K38_He4_to_Sc42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K38);
            }
#endif
            reaclib_rates::rate_K38_He4_to_Sc42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ca41_p_to_Sc42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca41);
            }
#endif
            reaclib_rates::rate_Ca41_p_to_Sc42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca41_He4_to_Ti45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca41);
            }
#endif
            reaclib_rates::rate_Ca41_He4_to_Ti45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc42_He4_to_V46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc42);
            }
#endif
            reaclib_rates::rate_Sc42_He4_to_V46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ti45_p_to_V46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti45);
            }
#endif
            reaclib_rates::rate_Ti45_p_to_V46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti45_He4_to_Cr49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti45);
            }
#endif
            reaclib_rates::rate_Ti45_He4_to_Cr49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V46_He4_to_Mn50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V46);
            }
#endif
            reaclib_rates::rate_V46_He4_to_Mn50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Cr49_p_to_Mn50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr49);
            }
#endif
            reaclib_rates::rate_Cr49_p_to_Mn50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // N13_He4_to_p_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N13);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N13);
            }
#endif
            reaclib_rates::rate_N13_He4_to_p_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N15_p_to_He4_C12_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_N15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N15);
            }
#endif
            reaclib_rates::rate_N15_p_to_He4_C12_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // O17_p_to_He4_N14_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O17);
            }
#endif
            reaclib_rates::rate_O17_p_to_He4_N14_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O18_p_to_He4_N15_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O18);
            }
#endif
            reaclib_rates::rate_O18_p_to_He4_N15_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F18_He4_to_p_Ne21_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F18);
            }
#endif
            reaclib_rates::rate_F18_He4_to_p_Ne21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F19_p_to_He4_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F19);
            }
#endif
            reaclib_rates::rate_F19_p_to_He4_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F19_He4_to_p_Ne22_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F19);
            }
#endif
            reaclib_rates::rate_F19_He4_to_p_Ne22_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Na22_He4_to_p_Mg25_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na22);
            }
#endif
            reaclib_rates::rate_Na22_He4_to_p_Mg25_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Al26_He4_to_p_Si29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al26);
            }
#endif
            reaclib_rates::rate_Al26_He4_to_p_Si29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // P30_He4_to_p_S33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P30);
            }
#endif
            reaclib_rates::rate_P30_He4_to_p_S33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Cl34_He4_to_p_Ar37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl34);
            }
#endif
            reaclib_rates::rate_Cl34_He4_to_p_Ar37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // K38_He4_to_p_Ca41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K38);
            }
#endif
            reaclib_rates::rate_K38_He4_to_p_Ca41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Sc42_He4_to_p_Ti45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc42);
            }
#endif
            reaclib_rates::rate_Sc42_He4_to_p_Ti45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // V46_He4_to_p_Cr49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V46);
            }
#endif
            reaclib_rates::rate_V46_He4_to_p_Cr49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Cr49_He4_to_p_Mn52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr49);
            }
#endif
            reaclib_rates::rate_Cr49_He4_to_p_Mn52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Mn50_n_to_Mn51_reaclib
            reaclib_rates::rate_Mn50_n_to_Mn51_reaclib(rate_eval, tfactors);
        }

        {
            // Mn51_n_to_Mn52_reaclib
            reaclib_rates::rate_Mn51_n_to_Mn52_reaclib(rate_eval, tfactors);
        }

        {
            // Mn52_p_to_Fe53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn52);
            }
#endif
            reaclib_rates::rate_Mn52_p_to_Fe53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn52_He4_to_Co56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn52);
            }
#endif
            reaclib_rates::rate_Mn52_He4_to_Co56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe52_n_to_Fe53_reaclib
            reaclib_rates::rate_Fe52_n_to_Fe53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_n_to_Fe54_reaclib
            reaclib_rates::rate_Fe53_n_to_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_He4_to_Ni57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe53);
            }
#endif
            reaclib_rates::rate_Fe53_He4_to_Ni57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe54_n_to_Fe55_reaclib
            reaclib_rates::rate_Fe54_n_to_Fe55_reaclib(rate_eval, tfactors);
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
            // Fe55_n_to_Fe56_reaclib
            reaclib_rates::rate_Fe55_n_to_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Fe55_p_to_Co56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe55);
            }
#endif
            reaclib_rates::rate_Fe55_p_to_Co56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Co55_n_to_Co56_reaclib
            reaclib_rates::rate_Co55_n_to_Co56_reaclib(rate_eval, tfactors);
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
            // Co56_p_to_Ni57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
            }
#endif
            reaclib_rates::rate_Co56_p_to_Ni57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ni56_n_to_Ni57_reaclib
            reaclib_rates::rate_Ni56_n_to_Ni57_reaclib(rate_eval, tfactors);
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
            // Ni57_n_to_Ni58_reaclib
            reaclib_rates::rate_Ni57_n_to_Ni58_reaclib(rate_eval, tfactors);
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
            // Mn50_He4_to_p_Fe53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn50);
            }
#endif
            reaclib_rates::rate_Mn50_He4_to_p_Fe53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Mn52_He4_to_p_Fe55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn52);
            }
#endif
            reaclib_rates::rate_Mn52_He4_to_p_Fe55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe52_n_to_p_Mn52_reaclib
            reaclib_rates::rate_Fe52_n_to_p_Mn52_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_He4_to_p_Co56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe53);
            }
#endif
            reaclib_rates::rate_Fe53_He4_to_p_Co56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_n_to_p_Fe55_reaclib
            reaclib_rates::rate_Co55_n_to_p_Fe55_reaclib(rate_eval, tfactors);
        }

        {
            // Co55_n_to_He4_Mn52_reaclib
            reaclib_rates::rate_Co55_n_to_He4_Mn52_reaclib(rate_eval, tfactors);
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
            // Ni56_n_to_p_Co56_reaclib
            reaclib_rates::rate_Ni56_n_to_p_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Ni56_n_to_He4_Fe53_reaclib
            reaclib_rates::rate_Ni56_n_to_He4_Fe53_reaclib(rate_eval, tfactors);
        }

        {
            // Ni57_n_to_p_Co57_reaclib
            reaclib_rates::rate_Ni57_n_to_p_Co57_reaclib(rate_eval, tfactors);
        }

        {
            // Ni57_n_to_He4_Fe54_reaclib
            reaclib_rates::rate_Ni57_n_to_He4_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Ni58_n_to_He4_Fe55_reaclib
            reaclib_rates::rate_Ni58_n_to_He4_Fe55_reaclib(rate_eval, tfactors);
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
            // Zn60_n_to_He4_Ni57_reaclib
            reaclib_rates::rate_Zn60_n_to_He4_Ni57_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_to_Mn52_electron_capture_weaktab
            tabular_weak_rates::rate_Fe52_to_Mn52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Fe52_beta_neg_weaktab
            tabular_weak_rates::rate_Mn52_to_Fe52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co55_to_Fe55_electron_capture_weaktab
            tabular_weak_rates::rate_Co55_to_Fe55_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe55_to_Co55_beta_neg_weaktab
            tabular_weak_rates::rate_Fe55_to_Co55_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Co57_to_Ni57_beta_neg_weaktab
            tabular_weak_rates::rate_Co57_to_Ni57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Co57_electron_capture_weaktab
            tabular_weak_rates::rate_Ni57_to_Co57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // F18_to_O18_electron_capture_weaktab
            tabular_weak_rates::rate_F18_to_O18_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // O18_to_F18_beta_neg_weaktab
            tabular_weak_rates::rate_O18_to_F18_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Na22_to_Ne22_electron_capture_weaktab
            tabular_weak_rates::rate_Na22_to_Ne22_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ne22_to_Na22_beta_neg_weaktab
            tabular_weak_rates::rate_Ne22_to_Na22_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // N13_to_p_C12_derived
            derived_rates::rate_N13_to_p_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N14_to_p_C13_derived
            derived_rates::rate_N14_to_p_C13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N15_to_p_C14_derived
            derived_rates::rate_N15_to_p_C14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O16_to_p_N15_derived
            derived_rates::rate_O16_to_p_N15_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O16_to_He4_C12_derived
            derived_rates::rate_O16_to_He4_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O18_to_He4_C14_derived
            derived_rates::rate_O18_to_He4_C14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // F18_to_p_O17_derived
            derived_rates::rate_F18_to_p_O17_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // F18_to_He4_N14_derived
            derived_rates::rate_F18_to_He4_N14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // F19_to_p_O18_derived
            derived_rates::rate_F19_to_p_O18_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // F19_to_He4_N15_derived
            derived_rates::rate_F19_to_He4_N15_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_to_p_F19_derived
            derived_rates::rate_Ne20_to_p_F19_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_to_He4_O16_derived
            derived_rates::rate_Ne20_to_He4_O16_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne21_to_He4_O17_derived
            derived_rates::rate_Ne21_to_He4_O17_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne22_to_He4_O18_derived
            derived_rates::rate_Ne22_to_He4_O18_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na22_to_p_Ne21_derived
            derived_rates::rate_Na22_to_p_Ne21_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na22_to_He4_F18_derived
            derived_rates::rate_Na22_to_He4_F18_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na23_to_p_Ne22_derived
            derived_rates::rate_Na23_to_p_Ne22_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na23_to_He4_F19_derived
            derived_rates::rate_Na23_to_He4_F19_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg24_to_p_Na23_derived
            derived_rates::rate_Mg24_to_p_Na23_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg24_to_He4_Ne20_derived
            derived_rates::rate_Mg24_to_He4_Ne20_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg25_to_He4_Ne21_derived
            derived_rates::rate_Mg25_to_He4_Ne21_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Al26_to_p_Mg25_derived
            derived_rates::rate_Al26_to_p_Mg25_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Al26_to_He4_Na22_derived
            derived_rates::rate_Al26_to_He4_Na22_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Al27_to_He4_Na23_derived
            derived_rates::rate_Al27_to_He4_Na23_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_to_p_Al27_derived
            derived_rates::rate_Si28_to_p_Al27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_to_He4_Mg24_derived
            derived_rates::rate_Si28_to_He4_Mg24_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si29_to_He4_Mg25_derived
            derived_rates::rate_Si29_to_He4_Mg25_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P30_to_p_Si29_derived
            derived_rates::rate_P30_to_p_Si29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P30_to_He4_Al26_derived
            derived_rates::rate_P30_to_He4_Al26_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P31_to_He4_Al27_derived
            derived_rates::rate_P31_to_He4_Al27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_to_p_P31_derived
            derived_rates::rate_S32_to_p_P31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_to_He4_Si28_derived
            derived_rates::rate_S32_to_He4_Si28_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S33_to_He4_Si29_derived
            derived_rates::rate_S33_to_He4_Si29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl34_to_p_S33_derived
            derived_rates::rate_Cl34_to_p_S33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl34_to_He4_P30_derived
            derived_rates::rate_Cl34_to_He4_P30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl35_to_He4_P31_derived
            derived_rates::rate_Cl35_to_He4_P31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar36_to_p_Cl35_derived
            derived_rates::rate_Ar36_to_p_Cl35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar36_to_He4_S32_derived
            derived_rates::rate_Ar36_to_He4_S32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar37_to_He4_S33_derived
            derived_rates::rate_Ar37_to_He4_S33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K38_to_p_Ar37_derived
            derived_rates::rate_K38_to_p_Ar37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K38_to_He4_Cl34_derived
            derived_rates::rate_K38_to_He4_Cl34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K39_to_He4_Cl35_derived
            derived_rates::rate_K39_to_He4_Cl35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_to_p_K39_derived
            derived_rates::rate_Ca40_to_p_K39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_to_He4_Ar36_derived
            derived_rates::rate_Ca40_to_He4_Ar36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca41_to_He4_Ar37_derived
            derived_rates::rate_Ca41_to_He4_Ar37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc42_to_p_Ca41_derived
            derived_rates::rate_Sc42_to_p_Ca41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc42_to_He4_K38_derived
            derived_rates::rate_Sc42_to_He4_K38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc43_to_He4_K39_derived
            derived_rates::rate_Sc43_to_He4_K39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti44_to_p_Sc43_derived
            derived_rates::rate_Ti44_to_p_Sc43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti44_to_He4_Ca40_derived
            derived_rates::rate_Ti44_to_He4_Ca40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti45_to_He4_Ca41_derived
            derived_rates::rate_Ti45_to_He4_Ca41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V46_to_p_Ti45_derived
            derived_rates::rate_V46_to_p_Ti45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V46_to_He4_Sc42_derived
            derived_rates::rate_V46_to_He4_Sc42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V47_to_He4_Sc43_derived
            derived_rates::rate_V47_to_He4_Sc43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr48_to_p_V47_derived
            derived_rates::rate_Cr48_to_p_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr48_to_He4_Ti44_derived
            derived_rates::rate_Cr48_to_He4_Ti44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr49_to_He4_Ti45_derived
            derived_rates::rate_Cr49_to_He4_Ti45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn50_to_p_Cr49_derived
            derived_rates::rate_Mn50_to_p_Cr49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn50_to_He4_V46_derived
            derived_rates::rate_Mn50_to_He4_V46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn51_to_He4_V47_derived
            derived_rates::rate_Mn51_to_He4_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe52_to_p_Mn51_derived
            derived_rates::rate_Fe52_to_p_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe52_to_He4_Cr48_derived
            derived_rates::rate_Fe52_to_He4_Cr48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_He4_Mn51_derived
            derived_rates::rate_Co55_to_He4_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni56_to_p_Co55_derived
            derived_rates::rate_Ni56_to_p_Co55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni56_to_He4_Fe52_derived
            derived_rates::rate_Ni56_to_He4_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C12_to_3He4_derived
            derived_rates::rate_C12_to_3He4_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C12_He4_to_p_N15_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C12);
            }
#endif
            derived_rates::rate_C12_He4_to_p_N15_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // N14_He4_to_p_O17_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N14);
            }
#endif
            derived_rates::rate_N14_He4_to_p_O17_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // N15_He4_to_p_O18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N15);
            }
#endif
            derived_rates::rate_N15_He4_to_p_O18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // O16_p_to_He4_N13_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O16);
            }
#endif
            derived_rates::rate_O16_p_to_He4_N13_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // O16_He4_to_p_F19_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O16);
            }
#endif
            derived_rates::rate_O16_He4_to_p_F19_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ne20_He4_to_p_Na23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            derived_rates::rate_Ne20_He4_to_p_Na23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ne20_He4_to_C12_C12_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            derived_rates::rate_Ne20_He4_to_C12_C12_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ne21_p_to_He4_F18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne21);
            }
#endif
            derived_rates::rate_Ne21_p_to_He4_F18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ne22_p_to_He4_F19_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne22);
            }
#endif
            derived_rates::rate_Ne22_p_to_He4_F19_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Na23_p_to_C12_C12_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            derived_rates::rate_Na23_p_to_C12_C12_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg24_He4_to_p_Al27_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            derived_rates::rate_Mg24_He4_to_p_Al27_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg24_He4_to_C12_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            derived_rates::rate_Mg24_He4_to_C12_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg25_p_to_He4_Na22_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg25);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg25);
            }
#endif
            derived_rates::rate_Mg25_p_to_He4_Na22_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Al27_p_to_C12_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            derived_rates::rate_Al27_p_to_C12_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si28_He4_to_p_P31_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            derived_rates::rate_Si28_He4_to_p_P31_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si28_He4_to_C12_Ne20_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            derived_rates::rate_Si28_He4_to_C12_Ne20_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si28_He4_to_O16_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            derived_rates::rate_Si28_He4_to_O16_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si29_p_to_He4_Al26_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si29);
            }
#endif
            derived_rates::rate_Si29_p_to_He4_Al26_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P31_p_to_C12_Ne20_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            derived_rates::rate_P31_p_to_C12_Ne20_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P31_p_to_O16_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            derived_rates::rate_P31_p_to_O16_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S32_He4_to_p_Cl35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
            }
#endif
            derived_rates::rate_S32_He4_to_p_Cl35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S33_p_to_He4_P30_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S33);
            }
#endif
            derived_rates::rate_S33_p_to_He4_P30_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar36_He4_to_p_K39_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
            }
#endif
            derived_rates::rate_Ar36_He4_to_p_K39_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar37_p_to_He4_Cl34_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar37);
            }
#endif
            derived_rates::rate_Ar37_p_to_He4_Cl34_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca40_He4_to_p_Sc43_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
            }
#endif
            derived_rates::rate_Ca40_He4_to_p_Sc43_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca41_p_to_He4_K38_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca41);
            }
#endif
            derived_rates::rate_Ca41_p_to_He4_K38_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti45_p_to_He4_Sc42_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti45);
            }
#endif
            derived_rates::rate_Ti45_p_to_He4_Sc42_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V47_p_to_He4_Ti44_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V47);
            }
#endif
            derived_rates::rate_V47_p_to_He4_Ti44_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr49_p_to_He4_V46_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr49);
            }
#endif
            derived_rates::rate_Cr49_p_to_He4_V46_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn51_p_to_He4_Cr48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
            }
#endif
            derived_rates::rate_Mn51_p_to_He4_Cr48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn52_p_to_He4_Cr49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn52);
            }
#endif
            derived_rates::rate_Mn52_p_to_He4_Cr49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co55_p_to_He4_Fe52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
            }
#endif
            derived_rates::rate_Co55_p_to_He4_Fe52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn51_to_n_Mn50_derived
            derived_rates::rate_Mn51_to_n_Mn50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn52_to_n_Mn51_derived
            derived_rates::rate_Mn52_to_n_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe53_to_n_Fe52_derived
            derived_rates::rate_Fe53_to_n_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe53_to_p_Mn52_derived
            derived_rates::rate_Fe53_to_p_Mn52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe54_to_n_Fe53_derived
            derived_rates::rate_Fe54_to_n_Fe53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe55_to_n_Fe54_derived
            derived_rates::rate_Fe55_to_n_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe56_to_n_Fe55_derived
            derived_rates::rate_Fe56_to_n_Fe55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_p_Fe54_derived
            derived_rates::rate_Co55_to_p_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co56_to_n_Co55_derived
            derived_rates::rate_Co56_to_n_Co55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co56_to_p_Fe55_derived
            derived_rates::rate_Co56_to_p_Fe55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co56_to_He4_Mn52_derived
            derived_rates::rate_Co56_to_He4_Mn52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co57_to_n_Co56_derived
            derived_rates::rate_Co57_to_n_Co56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co57_to_p_Fe56_derived
            derived_rates::rate_Co57_to_p_Fe56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni57_to_n_Ni56_derived
            derived_rates::rate_Ni57_to_n_Ni56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni57_to_p_Co56_derived
            derived_rates::rate_Ni57_to_p_Co56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni57_to_He4_Fe53_derived
            derived_rates::rate_Ni57_to_He4_Fe53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni58_to_n_Ni57_derived
            derived_rates::rate_Ni58_to_n_Ni57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni58_to_p_Co57_derived
            derived_rates::rate_Ni58_to_p_Co57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni58_to_He4_Fe54_derived
            derived_rates::rate_Ni58_to_He4_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu59_to_p_Ni58_derived
            derived_rates::rate_Cu59_to_p_Ni58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu59_to_He4_Co55_derived
            derived_rates::rate_Cu59_to_He4_Co55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn60_to_p_Cu59_derived
            derived_rates::rate_Zn60_to_p_Cu59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn60_to_He4_Ni56_derived
            derived_rates::rate_Zn60_to_He4_Ni56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn52_p_to_n_Fe52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn52);
            }
#endif
            derived_rates::rate_Mn52_p_to_n_Fe52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn52_He4_to_n_Co55_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn52);
            }
#endif
            derived_rates::rate_Mn52_He4_to_n_Co55_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe53_p_to_He4_Mn50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe53);
            }
#endif
            derived_rates::rate_Fe53_p_to_He4_Mn50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe53_He4_to_n_Ni56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe53);
            }
#endif
            derived_rates::rate_Fe53_He4_to_n_Ni56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe54_p_to_He4_Mn51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
            }
#endif
            derived_rates::rate_Fe54_p_to_He4_Mn51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe54_He4_to_n_Ni57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
            }
#endif
            derived_rates::rate_Fe54_He4_to_n_Ni57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe54_He4_to_p_Co57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
            }
#endif
            derived_rates::rate_Fe54_He4_to_p_Co57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe55_p_to_n_Co55_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe55);
            }
#endif
            derived_rates::rate_Fe55_p_to_n_Co55_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe55_p_to_He4_Mn52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe55);
            }
#endif
            derived_rates::rate_Fe55_p_to_He4_Mn52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe55_He4_to_n_Ni58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe55);
            }
#endif
            derived_rates::rate_Fe55_He4_to_n_Ni58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe56_p_to_n_Co56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            derived_rates::rate_Fe56_p_to_n_Co56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co56_p_to_n_Ni56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
            }
#endif
            derived_rates::rate_Co56_p_to_n_Ni56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co56_p_to_He4_Fe53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
            }
#endif
            derived_rates::rate_Co56_p_to_He4_Fe53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co56_He4_to_n_Cu59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co56);
            }
#endif
            derived_rates::rate_Co56_He4_to_n_Cu59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co57_p_to_n_Ni57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
            }
#endif
            derived_rates::rate_Co57_p_to_n_Ni57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni56_He4_to_p_Cu59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
            }
#endif
            derived_rates::rate_Ni56_He4_to_p_Cu59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni57_He4_to_n_Zn60_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni57);
            }
#endif
            derived_rates::rate_Ni57_He4_to_n_Zn60_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni58_p_to_He4_Co55_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
            }
#endif
            derived_rates::rate_Ni58_p_to_He4_Co55_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // C14_to_N14_beta_neg_reaclib
            reaclib_rates::rate_C14_to_N14_beta_neg_reaclib(rate_eval, tfactors);
        }

        {
            // N13_to_C13_beta_pos_reaclib
            reaclib_rates::rate_N13_to_C13_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_to_Mn52_electron_capture_weaktab
            tabular_weak_rates::rate_Fe52_to_Mn52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Fe52_beta_neg_weaktab
            tabular_weak_rates::rate_Mn52_to_Fe52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co55_to_Fe55_electron_capture_weaktab
            tabular_weak_rates::rate_Co55_to_Fe55_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe55_to_Co55_beta_neg_weaktab
            tabular_weak_rates::rate_Fe55_to_Co55_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Co57_to_Ni57_beta_neg_weaktab
            tabular_weak_rates::rate_Co57_to_Ni57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Co57_electron_capture_weaktab
            tabular_weak_rates::rate_Ni57_to_Co57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // F18_to_O18_electron_capture_weaktab
            tabular_weak_rates::rate_F18_to_O18_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // O18_to_F18_beta_neg_weaktab
            tabular_weak_rates::rate_O18_to_F18_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Na22_to_Ne22_electron_capture_weaktab
            tabular_weak_rates::rate_Na22_to_Ne22_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ne22_to_Na22_beta_neg_weaktab
            tabular_weak_rates::rate_Ne22_to_Na22_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }


    }

}

#endif
