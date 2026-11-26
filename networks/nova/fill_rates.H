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
            // N13_to_C13_beta_pos_reaclib
            reaclib_rates::rate_N13_to_C13_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // O14_to_N14_beta_pos_reaclib
            reaclib_rates::rate_O14_to_N14_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // O15_to_N15_beta_pos_reaclib
            reaclib_rates::rate_O15_to_N15_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // B8_to_He4_He4_beta_pos_reaclib
            reaclib_rates::rate_B8_to_He4_He4_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // p_p_to_d_beta_pos_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_p);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_p);
            }
#endif
            reaclib_rates::rate_p_p_to_d_beta_pos_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // p_p_to_d_electron_capture_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_p);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_p);
            }
#endif
            reaclib_rates::rate_p_p_to_d_electron_capture_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // d_p_to_He3_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_d);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_d);
            }
#endif
            reaclib_rates::rate_d_p_to_He3_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // d_d_to_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_d);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_d);
            }
#endif
            reaclib_rates::rate_d_d_to_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He3_p_to_He4_beta_pos_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_He3);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_He3);
            }
#endif
            reaclib_rates::rate_He3_p_to_He4_beta_pos_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He3_He4_to_Be7_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He3_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He3_He4);
            }
#endif
            reaclib_rates::rate_He3_He4_to_Be7_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Be7_p_to_B8_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Be7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Be7);
            }
#endif
            reaclib_rates::rate_Be7_p_to_B8_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // N13_p_to_O14_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_N13);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N13);
            }
#endif
            reaclib_rates::rate_N13_p_to_O14_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N14_p_to_O15_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_N14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N14);
            }
#endif
            reaclib_rates::rate_N14_p_to_O15_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // O16_p_to_F17_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O16);
            }
#endif
            reaclib_rates::rate_O16_p_to_F17_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // He3_d_to_p_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_He3);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_He3);
            }
#endif
            reaclib_rates::rate_He3_d_to_p_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // O14_He4_to_p_F17_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O14);
            }
#endif
            reaclib_rates::rate_O14_He4_to_p_F17_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // F18_p_to_He4_O15_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F18);
            }
#endif
            reaclib_rates::rate_F18_p_to_He4_O15_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He3_He3_to_p_p_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He3_He3);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He3_He3);
            }
#endif
            reaclib_rates::rate_He3_He3_to_p_p_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Be7_d_to_p_He4_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_Be7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_Be7);
            }
#endif
            reaclib_rates::rate_Be7_d_to_p_He4_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Be7_He3_to_p_p_He4_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He3_Be7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He3_Be7);
            }
#endif
            reaclib_rates::rate_Be7_He3_to_p_p_He4_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // F17_to_O17_electron_capture_weaktab
            tabular_weak_rates::rate_F17_to_O17_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // O17_to_F17_beta_neg_weaktab
            tabular_weak_rates::rate_O17_to_F17_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // He3_to_p_d_derived
            derived_rates::rate_He3_to_p_d_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // He4_to_d_d_derived
            derived_rates::rate_He4_to_d_d_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Be7_to_He4_He3_derived
            derived_rates::rate_Be7_to_He4_He3_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // B8_to_p_Be7_derived
            derived_rates::rate_B8_to_p_Be7_derived(rate_eval, tfactors, pf_cache);
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
            // O14_to_p_N13_derived
            derived_rates::rate_O14_to_p_N13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O15_to_p_N14_derived
            derived_rates::rate_O15_to_p_N14_derived(rate_eval, tfactors, pf_cache);
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
            // F17_to_p_O16_derived
            derived_rates::rate_F17_to_p_O16_derived(rate_eval, tfactors, pf_cache);
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
            // C12_to_3He4_derived
            derived_rates::rate_C12_to_3He4_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // He4_p_to_d_He3_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_He4);
            }
#endif
            derived_rates::rate_He4_p_to_d_He3_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // O15_He4_to_p_F18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O15);
            }
#endif
            derived_rates::rate_O15_He4_to_p_F18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // F17_p_to_He4_O14_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F17);
            }
#endif
            derived_rates::rate_F17_p_to_He4_O14_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_p_p_to_He3_He3_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_p) + rate_eval.log_screen(k_He2_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_p) + rate_eval.dlog_screen_dT(k_He2_He4);
            }
#endif
            derived_rates::rate_He4_p_p_to_He3_He3_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_He4_p_to_d_Be7_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_p_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_p_Be8);
            }
#endif
            derived_rates::rate_He4_He4_p_to_d_Be7_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_He4_p_p_to_He3_Be7_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_p_p) + rate_eval.log_screen(k_He2_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_p_p) + rate_eval.dlog_screen_dT(k_He2_Be8);
            }
#endif
            derived_rates::rate_He4_He4_p_p_to_He3_Be7_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // N13_to_C13_beta_pos_reaclib
            reaclib_rates::rate_N13_to_C13_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // O14_to_N14_beta_pos_reaclib
            reaclib_rates::rate_O14_to_N14_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // O15_to_N15_beta_pos_reaclib
            reaclib_rates::rate_O15_to_N15_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // B8_to_He4_He4_beta_pos_reaclib
            reaclib_rates::rate_B8_to_He4_He4_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // p_p_to_d_beta_pos_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_p);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_p);
            }
#endif
            reaclib_rates::rate_p_p_to_d_beta_pos_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // p_p_to_d_electron_capture_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_p);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_p);
            }
#endif
            reaclib_rates::rate_p_p_to_d_electron_capture_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He3_p_to_He4_beta_pos_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_He3);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_He3);
            }
#endif
            reaclib_rates::rate_He3_p_to_He4_beta_pos_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F17_to_O17_electron_capture_weaktab
            tabular_weak_rates::rate_F17_to_O17_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // O17_to_F17_beta_neg_weaktab
            tabular_weak_rates::rate_O17_to_F17_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }


    }

}

#endif
