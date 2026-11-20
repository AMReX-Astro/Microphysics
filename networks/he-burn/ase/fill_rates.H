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
            // O16_O16_to_n_S31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_O16_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_O16_O16);
            }
#endif
            reaclib_rates::rate_O16_O16_to_n_S31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ar36_to_He4_S32_derived
            derived_rates::rate_Ar36_to_He4_S32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar36_to_p_Cl35_derived
            derived_rates::rate_Ar36_to_p_Cl35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_He4_to_Ar36_approx
            approximate_rates::rate_S32_He4_to_Ar36_approx(rate_eval);
        }

        {
            // Ar36_to_S32_He4_approx
            approximate_rates::rate_Ar36_to_S32_He4_approx(rate_eval);
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
            // Ca40_to_He4_Ar36_derived
            derived_rates::rate_Ca40_to_He4_Ar36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_to_p_K39_derived
            derived_rates::rate_Ca40_to_p_K39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar36_He4_to_Ca40_approx
            approximate_rates::rate_Ar36_He4_to_Ca40_approx(rate_eval);
        }

        {
            // Ca40_to_Ar36_He4_approx
            approximate_rates::rate_Ca40_to_Ar36_He4_approx(rate_eval);
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
            // Ti44_to_He4_Ca40_derived
            derived_rates::rate_Ti44_to_He4_Ca40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti44_to_p_Sc43_derived
            derived_rates::rate_Ti44_to_p_Sc43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_He4_to_Ti44_approx
            approximate_rates::rate_Ca40_He4_to_Ti44_approx(rate_eval);
        }

        {
            // Ti44_to_Ca40_He4_approx
            approximate_rates::rate_Ti44_to_Ca40_He4_approx(rate_eval);
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
            // Cr48_to_He4_Ti44_derived
            derived_rates::rate_Cr48_to_He4_Ti44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr48_to_p_V47_derived
            derived_rates::rate_Cr48_to_p_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti44_He4_to_Cr48_approx
            approximate_rates::rate_Ti44_He4_to_Cr48_approx(rate_eval);
        }

        {
            // Cr48_to_Ti44_He4_approx
            approximate_rates::rate_Cr48_to_Ti44_He4_approx(rate_eval);
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
            // Fe52_to_He4_Cr48_derived
            derived_rates::rate_Fe52_to_He4_Cr48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe52_to_p_Mn51_derived
            derived_rates::rate_Fe52_to_p_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr48_He4_to_Fe52_approx
            approximate_rates::rate_Cr48_He4_to_Fe52_approx(rate_eval);
        }

        {
            // Fe52_to_Cr48_He4_approx
            approximate_rates::rate_Fe52_to_Cr48_He4_approx(rate_eval);
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
            // Ni56_to_He4_Fe52_derived
            derived_rates::rate_Ni56_to_He4_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni56_to_p_Co55_derived
            derived_rates::rate_Ni56_to_p_Co55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe52_He4_to_Ni56_approx
            approximate_rates::rate_Fe52_He4_to_Ni56_approx(rate_eval);
        }

        {
            // Ni56_to_Fe52_He4_approx
            approximate_rates::rate_Ni56_to_Fe52_He4_approx(rate_eval);
        }

        {
            // C12_C12_to_Mg24_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            modified_rates::rate_C12_C12_to_Mg24_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_O16_to_S32_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_O16_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_O16_O16);
            }
#endif
            modified_rates::rate_O16_O16_to_S32_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_C12_to_Si28_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            modified_rates::rate_O16_C12_to_Si28_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N13_to_p_C12_derived
            derived_rates::rate_N13_to_p_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O16_to_He4_C12_derived
            derived_rates::rate_O16_to_He4_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_to_He4_O16_derived
            derived_rates::rate_Ne20_to_He4_O16_derived(rate_eval, tfactors, pf_cache);
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
            // Si28_to_p_Al27_derived
            derived_rates::rate_Si28_to_p_Al27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_to_He4_Mg24_derived
            derived_rates::rate_Si28_to_He4_Mg24_derived(rate_eval, tfactors, pf_cache);
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
            // C12_to_3He4_derived
            derived_rates::rate_C12_to_3He4_derived(rate_eval, tfactors, pf_cache);
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
            // Mg24_to_C12_C12_derived
            derived_rates::rate_Mg24_to_C12_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_to_O16_O16_derived
            derived_rates::rate_S32_to_O16_O16_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_to_C12_O16_derived
            derived_rates::rate_Si28_to_C12_O16_derived(rate_eval, tfactors, pf_cache);
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


    }

}

#endif
