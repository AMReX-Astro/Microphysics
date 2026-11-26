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
            // Be7_to_Li7_electron_capture_reaclib
            reaclib_rates::rate_Be7_to_Li7_electron_capture_reaclib(rate_eval, tfactors);
        }

        {
            // C14_to_N14_beta_neg_reaclib
            reaclib_rates::rate_C14_to_N14_beta_neg_reaclib(rate_eval, tfactors);
        }

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
            // F17_to_O17_beta_pos_reaclib
            reaclib_rates::rate_F17_to_O17_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // F18_to_O18_beta_pos_reaclib
            reaclib_rates::rate_F18_to_O18_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Ne18_to_F18_beta_pos_reaclib
            reaclib_rates::rate_Ne18_to_F18_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Ne19_to_F19_beta_pos_reaclib
            reaclib_rates::rate_Ne19_to_F19_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Zn59_to_Cu59_beta_pos_reaclib
            reaclib_rates::rate_Zn59_to_Cu59_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Ge63_to_Ga63_beta_pos_reaclib
            reaclib_rates::rate_Ge63_to_Ga63_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // B8_to_He4_He4_beta_pos_reaclib
            reaclib_rates::rate_B8_to_He4_He4_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Zn59_to_p_Ni58_beta_pos_reaclib
            reaclib_rates::rate_Zn59_to_p_Ni58_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // p_n_to_d_reaclib
            reaclib_rates::rate_p_n_to_d_reaclib(rate_eval, tfactors);
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
            // d_He4_to_Li6_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_He4);
            }
#endif
            reaclib_rates::rate_d_He4_to_Li6_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He3_n_to_He4_reaclib
            reaclib_rates::rate_He3_n_to_He4_reaclib(rate_eval, tfactors);
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
            // Li6_n_to_Li7_reaclib
            reaclib_rates::rate_Li6_n_to_Li7_reaclib(rate_eval, tfactors);
        }

        {
            // Li6_p_to_Be7_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Li6);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Li6);
            }
#endif
            reaclib_rates::rate_Li6_p_to_Be7_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Li6_He4_to_B10_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Li6);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Li6);
            }
#endif
            reaclib_rates::rate_Li6_He4_to_B10_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Li7_He4_to_B11_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Li7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Li7);
            }
#endif
            reaclib_rates::rate_Li7_He4_to_B11_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Be9_p_to_B10_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Be9);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Be9);
            }
#endif
            reaclib_rates::rate_Be9_p_to_B10_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // B10_n_to_B11_reaclib
            reaclib_rates::rate_B10_n_to_B11_reaclib(rate_eval, tfactors);
        }

        {
            // B11_p_to_C12_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_B11);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_B11);
            }
#endif
            reaclib_rates::rate_B11_p_to_C12_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C12_n_to_C13_reaclib
            reaclib_rates::rate_C12_n_to_C13_reaclib(rate_eval, tfactors);
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
            // C13_n_to_C14_reaclib
            reaclib_rates::rate_C13_n_to_C14_reaclib(rate_eval, tfactors);
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
            // N13_n_to_N14_reaclib
            reaclib_rates::rate_N13_n_to_N14_reaclib(rate_eval, tfactors);
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
            // N14_n_to_N15_reaclib
            reaclib_rates::rate_N14_n_to_N15_reaclib(rate_eval, tfactors);
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
            // O14_n_to_O15_reaclib
            reaclib_rates::rate_O14_n_to_O15_reaclib(rate_eval, tfactors);
        }

        {
            // O14_He4_to_Ne18_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O14);
            }
#endif
            reaclib_rates::rate_O14_He4_to_Ne18_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O15_n_to_O16_reaclib
            reaclib_rates::rate_O15_n_to_O16_reaclib(rate_eval, tfactors);
        }

        {
            // O15_He4_to_Ne19_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O15);
            }
#endif
            reaclib_rates::rate_O15_He4_to_Ne19_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_n_to_O17_reaclib
            reaclib_rates::rate_O16_n_to_O17_reaclib(rate_eval, tfactors);
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
            // O17_n_to_O18_reaclib
            reaclib_rates::rate_O17_n_to_O18_reaclib(rate_eval, tfactors);
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
            // F17_n_to_F18_reaclib
            reaclib_rates::rate_F17_n_to_F18_reaclib(rate_eval, tfactors);
        }

        {
            // F17_p_to_Ne18_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F17);
            }
#endif
            reaclib_rates::rate_F17_p_to_Ne18_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F17_He4_to_Na21_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F17);
            }
#endif
            reaclib_rates::rate_F17_He4_to_Na21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F18_n_to_F19_reaclib
            reaclib_rates::rate_F18_n_to_F19_reaclib(rate_eval, tfactors);
        }

        {
            // F18_p_to_Ne19_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F18);
            }
#endif
            reaclib_rates::rate_F18_p_to_Ne19_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ne18_n_to_Ne19_reaclib
            reaclib_rates::rate_Ne18_n_to_Ne19_reaclib(rate_eval, tfactors);
        }

        {
            // Ne19_n_to_Ne20_reaclib
            reaclib_rates::rate_Ne19_n_to_Ne20_reaclib(rate_eval, tfactors);
        }

        {
            // Ne19_He4_to_Mg23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne19);
            }
#endif
            reaclib_rates::rate_Ne19_He4_to_Mg23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_n_to_Ne21_reaclib
            reaclib_rates::rate_Ne20_n_to_Ne21_reaclib(rate_eval, tfactors);
        }

        {
            // Ne20_p_to_Na21_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_p_to_Na21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ne21_n_to_Ne22_reaclib
            reaclib_rates::rate_Ne21_n_to_Ne22_reaclib(rate_eval, tfactors);
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
            // Ne22_He4_to_Mg26_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne22);
            }
#endif
            reaclib_rates::rate_Ne22_He4_to_Mg26_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na21_n_to_Na22_reaclib
            reaclib_rates::rate_Na21_n_to_Na22_reaclib(rate_eval, tfactors);
        }

        {
            // Na21_He4_to_Al25_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na21);
            }
#endif
            reaclib_rates::rate_Na21_He4_to_Al25_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na22_n_to_Na23_reaclib
            reaclib_rates::rate_Na22_n_to_Na23_reaclib(rate_eval, tfactors);
        }

        {
            // Na22_p_to_Mg23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na22);
            }
#endif
            reaclib_rates::rate_Na22_p_to_Mg23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Mg23_n_to_Mg24_reaclib
            reaclib_rates::rate_Mg23_n_to_Mg24_reaclib(rate_eval, tfactors);
        }

        {
            // Mg24_n_to_Mg25_reaclib
            reaclib_rates::rate_Mg24_n_to_Mg25_reaclib(rate_eval, tfactors);
        }

        {
            // Mg24_p_to_Al25_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_p_to_Al25_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Mg25_n_to_Mg26_reaclib
            reaclib_rates::rate_Mg25_n_to_Mg26_reaclib(rate_eval, tfactors);
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
            // Mg26_p_to_Al27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg26);
            }
#endif
            reaclib_rates::rate_Mg26_p_to_Al27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg26_He4_to_Si30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg26);
            }
#endif
            reaclib_rates::rate_Mg26_He4_to_Si30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al25_n_to_Al26_reaclib
            reaclib_rates::rate_Al25_n_to_Al26_reaclib(rate_eval, tfactors);
        }

        {
            // Al25_He4_to_P29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al25);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al25);
            }
#endif
            reaclib_rates::rate_Al25_He4_to_P29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al26_n_to_Al27_reaclib
            reaclib_rates::rate_Al26_n_to_Al27_reaclib(rate_eval, tfactors);
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
            // Si31_n_to_Si32_reaclib
            reaclib_rates::rate_Si31_n_to_Si32_reaclib(rate_eval, tfactors);
        }

        {
            // Si31_p_to_P32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si31);
            }
#endif
            reaclib_rates::rate_Si31_p_to_P32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Si32_p_to_P33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si32);
            }
#endif
            reaclib_rates::rate_Si32_p_to_P33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si32_He4_to_S36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si32);
            }
#endif
            reaclib_rates::rate_Si32_He4_to_S36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P29_n_to_P30_reaclib
            reaclib_rates::rate_P29_n_to_P30_reaclib(rate_eval, tfactors);
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
            // P30_n_to_P31_reaclib
            reaclib_rates::rate_P30_n_to_P31_reaclib(rate_eval, tfactors);
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
            // P31_n_to_P32_reaclib
            reaclib_rates::rate_P31_n_to_P32_reaclib(rate_eval, tfactors);
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
            // P32_n_to_P33_reaclib
            reaclib_rates::rate_P32_n_to_P33_reaclib(rate_eval, tfactors);
        }

        {
            // P32_p_to_S33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P32);
            }
#endif
            reaclib_rates::rate_P32_p_to_S33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P32_He4_to_Cl36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P32);
            }
#endif
            reaclib_rates::rate_P32_He4_to_Cl36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // S35_n_to_S36_reaclib
            reaclib_rates::rate_S35_n_to_S36_reaclib(rate_eval, tfactors);
        }

        {
            // S35_p_to_Cl36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S35);
            }
#endif
            reaclib_rates::rate_S35_p_to_Cl36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S35_He4_to_Ar39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S35);
            }
#endif
            reaclib_rates::rate_S35_He4_to_Ar39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S36_p_to_Cl37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S36);
            }
#endif
            reaclib_rates::rate_S36_p_to_Cl37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S36_He4_to_Ar40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S36);
            }
#endif
            reaclib_rates::rate_S36_He4_to_Ar40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl33_n_to_Cl34_reaclib
            reaclib_rates::rate_Cl33_n_to_Cl34_reaclib(rate_eval, tfactors);
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
            // Cl34_n_to_Cl35_reaclib
            reaclib_rates::rate_Cl34_n_to_Cl35_reaclib(rate_eval, tfactors);
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
            // Cl35_n_to_Cl36_reaclib
            reaclib_rates::rate_Cl35_n_to_Cl36_reaclib(rate_eval, tfactors);
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
            // Cl36_n_to_Cl37_reaclib
            reaclib_rates::rate_Cl36_n_to_Cl37_reaclib(rate_eval, tfactors);
        }

        {
            // Cl36_p_to_Ar37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl36);
            }
#endif
            reaclib_rates::rate_Cl36_p_to_Ar37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl36_He4_to_K40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl36);
            }
#endif
            reaclib_rates::rate_Cl36_He4_to_K40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Cl37_He4_to_K41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl37);
            }
#endif
            reaclib_rates::rate_Cl37_He4_to_K41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ar38_n_to_Ar39_reaclib
            reaclib_rates::rate_Ar38_n_to_Ar39_reaclib(rate_eval, tfactors);
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
            // Ar39_n_to_Ar40_reaclib
            reaclib_rates::rate_Ar39_n_to_Ar40_reaclib(rate_eval, tfactors);
        }

        {
            // Ar39_p_to_K40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar39);
            }
#endif
            reaclib_rates::rate_Ar39_p_to_K40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar39_He4_to_Ca43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar39);
            }
#endif
            reaclib_rates::rate_Ar39_He4_to_Ca43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar40_p_to_K41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar40);
            }
#endif
            reaclib_rates::rate_Ar40_p_to_K41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar40_He4_to_Ca44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar40);
            }
#endif
            reaclib_rates::rate_Ar40_He4_to_Ca44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K37_n_to_K38_reaclib
            reaclib_rates::rate_K37_n_to_K38_reaclib(rate_eval, tfactors);
        }

        {
            // K38_n_to_K39_reaclib
            reaclib_rates::rate_K38_n_to_K39_reaclib(rate_eval, tfactors);
        }

        {
            // K39_n_to_K40_reaclib
            reaclib_rates::rate_K39_n_to_K40_reaclib(rate_eval, tfactors);
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
            // K40_n_to_K41_reaclib
            reaclib_rates::rate_K40_n_to_K41_reaclib(rate_eval, tfactors);
        }

        {
            // K40_p_to_Ca41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K40);
            }
#endif
            reaclib_rates::rate_K40_p_to_Ca41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K40_He4_to_Sc44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K40);
            }
#endif
            reaclib_rates::rate_K40_He4_to_Sc44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K41_p_to_Ca42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K41);
            }
#endif
            reaclib_rates::rate_K41_p_to_Ca42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K41_He4_to_Sc45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K41);
            }
#endif
            reaclib_rates::rate_K41_He4_to_Sc45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_n_to_Ca41_reaclib
            reaclib_rates::rate_Ca40_n_to_Ca41_reaclib(rate_eval, tfactors);
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
            // Ca41_n_to_Ca42_reaclib
            reaclib_rates::rate_Ca41_n_to_Ca42_reaclib(rate_eval, tfactors);
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
            // Ca42_n_to_Ca43_reaclib
            reaclib_rates::rate_Ca42_n_to_Ca43_reaclib(rate_eval, tfactors);
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
            // Ca43_n_to_Ca44_reaclib
            reaclib_rates::rate_Ca43_n_to_Ca44_reaclib(rate_eval, tfactors);
        }

        {
            // Ca43_p_to_Sc44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca43);
            }
#endif
            reaclib_rates::rate_Ca43_p_to_Sc44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca43_He4_to_Ti47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca43);
            }
#endif
            reaclib_rates::rate_Ca43_He4_to_Ti47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca44_n_to_Ca45_reaclib
            reaclib_rates::rate_Ca44_n_to_Ca45_reaclib(rate_eval, tfactors);
        }

        {
            // Ca44_p_to_Sc45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca44);
            }
#endif
            reaclib_rates::rate_Ca44_p_to_Sc45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca44_He4_to_Ti48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca44);
            }
#endif
            reaclib_rates::rate_Ca44_He4_to_Ti48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca45_n_to_Ca46_reaclib
            reaclib_rates::rate_Ca45_n_to_Ca46_reaclib(rate_eval, tfactors);
        }

        {
            // Ca45_p_to_Sc46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca45);
            }
#endif
            reaclib_rates::rate_Ca45_p_to_Sc46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca45_He4_to_Ti49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca45);
            }
#endif
            reaclib_rates::rate_Ca45_He4_to_Ti49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca46_n_to_Ca47_reaclib
            reaclib_rates::rate_Ca46_n_to_Ca47_reaclib(rate_eval, tfactors);
        }

        {
            // Ca46_p_to_Sc47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca46);
            }
#endif
            reaclib_rates::rate_Ca46_p_to_Sc47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca46_He4_to_Ti50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca46);
            }
#endif
            reaclib_rates::rate_Ca46_He4_to_Ti50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca47_n_to_Ca48_reaclib
            reaclib_rates::rate_Ca47_n_to_Ca48_reaclib(rate_eval, tfactors);
        }

        {
            // Ca47_p_to_Sc48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca47);
            }
#endif
            reaclib_rates::rate_Ca47_p_to_Sc48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca47_He4_to_Ti51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca47);
            }
#endif
            reaclib_rates::rate_Ca47_He4_to_Ti51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca48_p_to_Sc49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca48);
            }
#endif
            reaclib_rates::rate_Ca48_p_to_Sc49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_n_to_Sc44_reaclib
            reaclib_rates::rate_Sc43_n_to_Sc44_reaclib(rate_eval, tfactors);
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
            // Sc44_n_to_Sc45_reaclib
            reaclib_rates::rate_Sc44_n_to_Sc45_reaclib(rate_eval, tfactors);
        }

        {
            // Sc44_p_to_Ti45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc44);
            }
#endif
            reaclib_rates::rate_Sc44_p_to_Ti45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc44_He4_to_V48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc44);
            }
#endif
            reaclib_rates::rate_Sc44_He4_to_V48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc45_n_to_Sc46_reaclib
            reaclib_rates::rate_Sc45_n_to_Sc46_reaclib(rate_eval, tfactors);
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
            // Sc46_n_to_Sc47_reaclib
            reaclib_rates::rate_Sc46_n_to_Sc47_reaclib(rate_eval, tfactors);
        }

        {
            // Sc46_p_to_Ti47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc46);
            }
#endif
            reaclib_rates::rate_Sc46_p_to_Ti47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc46_He4_to_V50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc46);
            }
#endif
            reaclib_rates::rate_Sc46_He4_to_V50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc47_n_to_Sc48_reaclib
            reaclib_rates::rate_Sc47_n_to_Sc48_reaclib(rate_eval, tfactors);
        }

        {
            // Sc47_p_to_Ti48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc47);
            }
#endif
            reaclib_rates::rate_Sc47_p_to_Ti48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc47_He4_to_V51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc47);
            }
#endif
            reaclib_rates::rate_Sc47_He4_to_V51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc48_n_to_Sc49_reaclib
            reaclib_rates::rate_Sc48_n_to_Sc49_reaclib(rate_eval, tfactors);
        }

        {
            // Sc48_p_to_Ti49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc48);
            }
#endif
            reaclib_rates::rate_Sc48_p_to_Ti49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc48_He4_to_V52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc48);
            }
#endif
            reaclib_rates::rate_Sc48_He4_to_V52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc49_p_to_Ti50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc49);
            }
#endif
            reaclib_rates::rate_Sc49_p_to_Ti50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_n_to_Ti45_reaclib
            reaclib_rates::rate_Ti44_n_to_Ti45_reaclib(rate_eval, tfactors);
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
            // Ti45_n_to_Ti46_reaclib
            reaclib_rates::rate_Ti45_n_to_Ti46_reaclib(rate_eval, tfactors);
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
            // Ti46_n_to_Ti47_reaclib
            reaclib_rates::rate_Ti46_n_to_Ti47_reaclib(rate_eval, tfactors);
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
            // Ti47_n_to_Ti48_reaclib
            reaclib_rates::rate_Ti47_n_to_Ti48_reaclib(rate_eval, tfactors);
        }

        {
            // Ti47_p_to_V48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti47);
            }
#endif
            reaclib_rates::rate_Ti47_p_to_V48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti47_He4_to_Cr51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti47);
            }
#endif
            reaclib_rates::rate_Ti47_He4_to_Cr51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti48_n_to_Ti49_reaclib
            reaclib_rates::rate_Ti48_n_to_Ti49_reaclib(rate_eval, tfactors);
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
            // Ti49_n_to_Ti50_reaclib
            reaclib_rates::rate_Ti49_n_to_Ti50_reaclib(rate_eval, tfactors);
        }

        {
            // Ti49_p_to_V50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti49);
            }
#endif
            reaclib_rates::rate_Ti49_p_to_V50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti49_He4_to_Cr53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti49);
            }
#endif
            reaclib_rates::rate_Ti49_He4_to_Cr53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti50_n_to_Ti51_reaclib
            reaclib_rates::rate_Ti50_n_to_Ti51_reaclib(rate_eval, tfactors);
        }

        {
            // Ti50_p_to_V51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti50);
            }
#endif
            reaclib_rates::rate_Ti50_p_to_V51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti50_He4_to_Cr54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti50);
            }
#endif
            reaclib_rates::rate_Ti50_He4_to_Cr54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti51_p_to_V52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti51);
            }
#endif
            reaclib_rates::rate_Ti51_p_to_V52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V46_n_to_V47_reaclib
            reaclib_rates::rate_V46_n_to_V47_reaclib(rate_eval, tfactors);
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
            // V47_n_to_V48_reaclib
            reaclib_rates::rate_V47_n_to_V48_reaclib(rate_eval, tfactors);
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
            // V48_n_to_V49_reaclib
            reaclib_rates::rate_V48_n_to_V49_reaclib(rate_eval, tfactors);
        }

        {
            // V48_p_to_Cr49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V48);
            }
#endif
            reaclib_rates::rate_V48_p_to_Cr49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V48_He4_to_Mn52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V48);
            }
#endif
            reaclib_rates::rate_V48_He4_to_Mn52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_n_to_V50_reaclib
            reaclib_rates::rate_V49_n_to_V50_reaclib(rate_eval, tfactors);
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
            // V50_n_to_V51_reaclib
            reaclib_rates::rate_V50_n_to_V51_reaclib(rate_eval, tfactors);
        }

        {
            // V50_p_to_Cr51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V50);
            }
#endif
            reaclib_rates::rate_V50_p_to_Cr51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V50_He4_to_Mn54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V50);
            }
#endif
            reaclib_rates::rate_V50_He4_to_Mn54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V51_n_to_V52_reaclib
            reaclib_rates::rate_V51_n_to_V52_reaclib(rate_eval, tfactors);
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
            // V52_p_to_Cr53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V52);
            }
#endif
            reaclib_rates::rate_V52_p_to_Cr53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr48_n_to_Cr49_reaclib
            reaclib_rates::rate_Cr48_n_to_Cr49_reaclib(rate_eval, tfactors);
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
            // Cr49_n_to_Cr50_reaclib
            reaclib_rates::rate_Cr49_n_to_Cr50_reaclib(rate_eval, tfactors);
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
            // Cr49_He4_to_Fe53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr49);
            }
#endif
            reaclib_rates::rate_Cr49_He4_to_Fe53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr50_n_to_Cr51_reaclib
            reaclib_rates::rate_Cr50_n_to_Cr51_reaclib(rate_eval, tfactors);
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
            // Cr51_n_to_Cr52_reaclib
            reaclib_rates::rate_Cr51_n_to_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // Cr51_p_to_Mn52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr51);
            }
#endif
            reaclib_rates::rate_Cr51_p_to_Mn52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr51_He4_to_Fe55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr51);
            }
#endif
            reaclib_rates::rate_Cr51_He4_to_Fe55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr52_n_to_Cr53_reaclib
            reaclib_rates::rate_Cr52_n_to_Cr53_reaclib(rate_eval, tfactors);
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
            // Cr53_n_to_Cr54_reaclib
            reaclib_rates::rate_Cr53_n_to_Cr54_reaclib(rate_eval, tfactors);
        }

        {
            // Cr53_p_to_Mn54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr53);
            }
#endif
            reaclib_rates::rate_Cr53_p_to_Mn54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr53_He4_to_Fe57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr53);
            }
#endif
            reaclib_rates::rate_Cr53_He4_to_Fe57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Mn50_n_to_Mn51_reaclib
            reaclib_rates::rate_Mn50_n_to_Mn51_reaclib(rate_eval, tfactors);
        }

        {
            // Mn50_He4_to_Co54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn50);
            }
#endif
            reaclib_rates::rate_Mn50_He4_to_Co54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_n_to_Mn52_reaclib
            reaclib_rates::rate_Mn51_n_to_Mn52_reaclib(rate_eval, tfactors);
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
            // Mn52_n_to_Mn53_reaclib
            reaclib_rates::rate_Mn52_n_to_Mn53_reaclib(rate_eval, tfactors);
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
            // Mn53_n_to_Mn54_reaclib
            reaclib_rates::rate_Mn53_n_to_Mn54_reaclib(rate_eval, tfactors);
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
            // Mn54_n_to_Mn55_reaclib
            reaclib_rates::rate_Mn54_n_to_Mn55_reaclib(rate_eval, tfactors);
        }

        {
            // Mn54_p_to_Fe55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn54);
            }
#endif
            reaclib_rates::rate_Mn54_p_to_Fe55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn54_He4_to_Co58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn54);
            }
#endif
            reaclib_rates::rate_Mn54_He4_to_Co58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Fe52_n_to_Fe53_reaclib
            reaclib_rates::rate_Fe52_n_to_Fe53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_p_to_Co53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe52);
            }
#endif
            reaclib_rates::rate_Fe52_p_to_Co53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Fe53_n_to_Fe54_reaclib
            reaclib_rates::rate_Fe53_n_to_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_p_to_Co54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe53);
            }
#endif
            reaclib_rates::rate_Fe53_p_to_Co54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Fe55_He4_to_Ni59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe55);
            }
#endif
            reaclib_rates::rate_Fe55_He4_to_Ni59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe56_n_to_Fe57_reaclib
            reaclib_rates::rate_Fe56_n_to_Fe57_reaclib(rate_eval, tfactors);
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
            // Fe57_n_to_Fe58_reaclib
            reaclib_rates::rate_Fe57_n_to_Fe58_reaclib(rate_eval, tfactors);
        }

        {
            // Fe57_p_to_Co58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe57);
            }
#endif
            reaclib_rates::rate_Fe57_p_to_Co58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe57_He4_to_Ni61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe57);
            }
#endif
            reaclib_rates::rate_Fe57_He4_to_Ni61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Fe58_He4_to_Ni62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe58);
            }
#endif
            reaclib_rates::rate_Fe58_He4_to_Ni62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co53_n_to_Co54_reaclib
            reaclib_rates::rate_Co53_n_to_Co54_reaclib(rate_eval, tfactors);
        }

        {
            // Co53_He4_to_Cu57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co53);
            }
#endif
            reaclib_rates::rate_Co53_He4_to_Cu57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co54_n_to_Co55_reaclib
            reaclib_rates::rate_Co54_n_to_Co55_reaclib(rate_eval, tfactors);
        }

        {
            // Co54_He4_to_Cu58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co54);
            }
#endif
            reaclib_rates::rate_Co54_He4_to_Cu58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Co56_He4_to_Cu60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co56);
            }
#endif
            reaclib_rates::rate_Co56_He4_to_Cu60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co57_n_to_Co58_reaclib
            reaclib_rates::rate_Co57_n_to_Co58_reaclib(rate_eval, tfactors);
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
            // Co57_He4_to_Cu61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co57);
            }
#endif
            reaclib_rates::rate_Co57_He4_to_Cu61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co58_n_to_Co59_reaclib
            reaclib_rates::rate_Co58_n_to_Co59_reaclib(rate_eval, tfactors);
        }

        {
            // Co58_p_to_Ni59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co58);
            }
#endif
            reaclib_rates::rate_Co58_p_to_Ni59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co58_He4_to_Cu62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co58);
            }
#endif
            reaclib_rates::rate_Co58_He4_to_Cu62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Co59_He4_to_Cu63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co59);
            }
#endif
            reaclib_rates::rate_Co59_He4_to_Cu63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni56_n_to_Ni57_reaclib
            reaclib_rates::rate_Ni56_n_to_Ni57_reaclib(rate_eval, tfactors);
        }

        {
            // Ni56_p_to_Cu57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni56);
            }
#endif
            reaclib_rates::rate_Ni56_p_to_Cu57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ni57_p_to_Cu58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni57);
            }
#endif
            reaclib_rates::rate_Ni57_p_to_Cu58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni57_He4_to_Zn61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni57);
            }
#endif
            reaclib_rates::rate_Ni57_He4_to_Zn61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni58_n_to_Ni59_reaclib
            reaclib_rates::rate_Ni58_n_to_Ni59_reaclib(rate_eval, tfactors);
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
            // Ni58_He4_to_Zn62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni58);
            }
#endif
            reaclib_rates::rate_Ni58_He4_to_Zn62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni59_n_to_Ni60_reaclib
            reaclib_rates::rate_Ni59_n_to_Ni60_reaclib(rate_eval, tfactors);
        }

        {
            // Ni59_p_to_Cu60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni59);
            }
#endif
            reaclib_rates::rate_Ni59_p_to_Cu60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni59_He4_to_Zn63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni59);
            }
#endif
            reaclib_rates::rate_Ni59_He4_to_Zn63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni60_n_to_Ni61_reaclib
            reaclib_rates::rate_Ni60_n_to_Ni61_reaclib(rate_eval, tfactors);
        }

        {
            // Ni60_p_to_Cu61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni60);
            }
#endif
            reaclib_rates::rate_Ni60_p_to_Cu61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni60_He4_to_Zn64_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni60);
            }
#endif
            reaclib_rates::rate_Ni60_He4_to_Zn64_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni61_n_to_Ni62_reaclib
            reaclib_rates::rate_Ni61_n_to_Ni62_reaclib(rate_eval, tfactors);
        }

        {
            // Ni61_p_to_Cu62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni61);
            }
#endif
            reaclib_rates::rate_Ni61_p_to_Cu62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni61_He4_to_Zn65_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni61);
            }
#endif
            reaclib_rates::rate_Ni61_He4_to_Zn65_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni62_n_to_Ni63_reaclib
            reaclib_rates::rate_Ni62_n_to_Ni63_reaclib(rate_eval, tfactors);
        }

        {
            // Ni62_p_to_Cu63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni62);
            }
#endif
            reaclib_rates::rate_Ni62_p_to_Cu63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni62_He4_to_Zn66_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni62);
            }
#endif
            reaclib_rates::rate_Ni62_He4_to_Zn66_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni63_n_to_Ni64_reaclib
            reaclib_rates::rate_Ni63_n_to_Ni64_reaclib(rate_eval, tfactors);
        }

        {
            // Ni63_p_to_Cu64_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni63);
            }
#endif
            reaclib_rates::rate_Ni63_p_to_Cu64_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni64_p_to_Cu65_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni64);
            }
#endif
            reaclib_rates::rate_Ni64_p_to_Cu65_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu57_n_to_Cu58_reaclib
            reaclib_rates::rate_Cu57_n_to_Cu58_reaclib(rate_eval, tfactors);
        }

        {
            // Cu58_n_to_Cu59_reaclib
            reaclib_rates::rate_Cu58_n_to_Cu59_reaclib(rate_eval, tfactors);
        }

        {
            // Cu58_p_to_Zn59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu58);
            }
#endif
            reaclib_rates::rate_Cu58_p_to_Zn59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu58_He4_to_Ga62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu58);
            }
#endif
            reaclib_rates::rate_Cu58_He4_to_Ga62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu59_n_to_Cu60_reaclib
            reaclib_rates::rate_Cu59_n_to_Cu60_reaclib(rate_eval, tfactors);
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
            // Cu59_He4_to_Ga63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu59);
            }
#endif
            reaclib_rates::rate_Cu59_He4_to_Ga63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu60_n_to_Cu61_reaclib
            reaclib_rates::rate_Cu60_n_to_Cu61_reaclib(rate_eval, tfactors);
        }

        {
            // Cu60_p_to_Zn61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu60);
            }
#endif
            reaclib_rates::rate_Cu60_p_to_Zn61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu60_He4_to_Ga64_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu60);
            }
#endif
            reaclib_rates::rate_Cu60_He4_to_Ga64_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu61_n_to_Cu62_reaclib
            reaclib_rates::rate_Cu61_n_to_Cu62_reaclib(rate_eval, tfactors);
        }

        {
            // Cu61_p_to_Zn62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu61);
            }
#endif
            reaclib_rates::rate_Cu61_p_to_Zn62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu62_n_to_Cu63_reaclib
            reaclib_rates::rate_Cu62_n_to_Cu63_reaclib(rate_eval, tfactors);
        }

        {
            // Cu62_p_to_Zn63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu62);
            }
#endif
            reaclib_rates::rate_Cu62_p_to_Zn63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu63_n_to_Cu64_reaclib
            reaclib_rates::rate_Cu63_n_to_Cu64_reaclib(rate_eval, tfactors);
        }

        {
            // Cu63_p_to_Zn64_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu63);
            }
#endif
            reaclib_rates::rate_Cu63_p_to_Zn64_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu64_n_to_Cu65_reaclib
            reaclib_rates::rate_Cu64_n_to_Cu65_reaclib(rate_eval, tfactors);
        }

        {
            // Cu64_p_to_Zn65_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu64);
            }
#endif
            reaclib_rates::rate_Cu64_p_to_Zn65_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu65_p_to_Zn66_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu65);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu65);
            }
#endif
            reaclib_rates::rate_Cu65_p_to_Zn66_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn59_n_to_Zn60_reaclib
            reaclib_rates::rate_Zn59_n_to_Zn60_reaclib(rate_eval, tfactors);
        }

        {
            // Zn59_He4_to_Ge63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Zn59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Zn59);
            }
#endif
            reaclib_rates::rate_Zn59_He4_to_Ge63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn60_n_to_Zn61_reaclib
            reaclib_rates::rate_Zn60_n_to_Zn61_reaclib(rate_eval, tfactors);
        }

        {
            // Zn60_He4_to_Ge64_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Zn60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Zn60);
            }
#endif
            reaclib_rates::rate_Zn60_He4_to_Ge64_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn61_n_to_Zn62_reaclib
            reaclib_rates::rate_Zn61_n_to_Zn62_reaclib(rate_eval, tfactors);
        }

        {
            // Zn61_p_to_Ga62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn61);
            }
#endif
            reaclib_rates::rate_Zn61_p_to_Ga62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn62_n_to_Zn63_reaclib
            reaclib_rates::rate_Zn62_n_to_Zn63_reaclib(rate_eval, tfactors);
        }

        {
            // Zn62_p_to_Ga63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn62);
            }
#endif
            reaclib_rates::rate_Zn62_p_to_Ga63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn63_n_to_Zn64_reaclib
            reaclib_rates::rate_Zn63_n_to_Zn64_reaclib(rate_eval, tfactors);
        }

        {
            // Zn63_p_to_Ga64_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn63);
            }
#endif
            reaclib_rates::rate_Zn63_p_to_Ga64_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn64_n_to_Zn65_reaclib
            reaclib_rates::rate_Zn64_n_to_Zn65_reaclib(rate_eval, tfactors);
        }

        {
            // Zn65_n_to_Zn66_reaclib
            reaclib_rates::rate_Zn65_n_to_Zn66_reaclib(rate_eval, tfactors);
        }

        {
            // Ga62_n_to_Ga63_reaclib
            reaclib_rates::rate_Ga62_n_to_Ga63_reaclib(rate_eval, tfactors);
        }

        {
            // Ga62_p_to_Ge63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ga62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ga62);
            }
#endif
            reaclib_rates::rate_Ga62_p_to_Ge63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ga63_n_to_Ga64_reaclib
            reaclib_rates::rate_Ga63_n_to_Ga64_reaclib(rate_eval, tfactors);
        }

        {
            // Ga63_p_to_Ge64_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ga63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ga63);
            }
#endif
            reaclib_rates::rate_Ga63_p_to_Ge64_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ge63_n_to_Ge64_reaclib
            reaclib_rates::rate_Ge63_n_to_Ge64_reaclib(rate_eval, tfactors);
        }

        {
            // d_d_to_n_He3_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_d);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_d);
            }
#endif
            reaclib_rates::rate_d_d_to_n_He3_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Li6_p_to_He4_He3_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Li6);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Li6);
            }
#endif
            reaclib_rates::rate_Li6_p_to_He4_He3_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Li6_d_to_n_Be7_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_Li6);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_Li6);
            }
#endif
            reaclib_rates::rate_Li6_d_to_n_Be7_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Li6_d_to_p_Li7_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_Li6);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_Li6);
            }
#endif
            reaclib_rates::rate_Li6_d_to_p_Li7_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Li7_p_to_He4_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Li7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Li7);
            }
#endif
            reaclib_rates::rate_Li7_p_to_He4_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Be7_n_to_p_Li7_reaclib
            reaclib_rates::rate_Be7_n_to_p_Li7_reaclib(rate_eval, tfactors);
        }

        {
            // Be7_n_to_He4_He4_reaclib
            reaclib_rates::rate_Be7_n_to_He4_He4_reaclib(rate_eval, tfactors);
        }

        {
            // Be9_p_to_He4_Li6_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Be9);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Be9);
            }
#endif
            reaclib_rates::rate_Be9_p_to_He4_Li6_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Be9_He4_to_n_C12_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Be9);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Be9);
            }
#endif
            reaclib_rates::rate_Be9_He4_to_n_C12_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // B10_n_to_He4_Li7_reaclib
            reaclib_rates::rate_B10_n_to_He4_Li7_reaclib(rate_eval, tfactors);
        }

        {
            // B10_p_to_He4_Be7_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_B10);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_B10);
            }
#endif
            reaclib_rates::rate_B10_p_to_He4_Be7_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // B10_He4_to_n_N13_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_B10);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_B10);
            }
#endif
            reaclib_rates::rate_B10_He4_to_n_N13_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // B10_He4_to_p_C13_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_B10);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_B10);
            }
#endif
            reaclib_rates::rate_B10_He4_to_p_C13_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // B11_He4_to_n_N14_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_B11);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_B11);
            }
#endif
            reaclib_rates::rate_B11_He4_to_n_N14_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // B11_He4_to_p_C14_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_B11);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_B11);
            }
#endif
            reaclib_rates::rate_B11_He4_to_p_C14_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // C13_p_to_n_N13_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_C13);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C13);
            }
#endif
            reaclib_rates::rate_C13_p_to_n_N13_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C13_d_to_n_N14_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_C13);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_C13);
            }
#endif
            reaclib_rates::rate_C13_d_to_n_N14_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C13_He4_to_n_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C13);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C13);
            }
#endif
            reaclib_rates::rate_C13_He4_to_n_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // C14_d_to_n_N15_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_C14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_C14);
            }
#endif
            reaclib_rates::rate_C14_d_to_n_N15_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // N14_n_to_p_C14_reaclib
            reaclib_rates::rate_N14_n_to_p_C14_reaclib(rate_eval, tfactors);
        }

        {
            // N14_p_to_n_O14_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_N14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N14);
            }
#endif
            reaclib_rates::rate_N14_p_to_n_O14_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N14_He4_to_n_F17_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N14);
            }
#endif
            reaclib_rates::rate_N14_He4_to_n_F17_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // N15_p_to_n_O15_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_N15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N15);
            }
#endif
            reaclib_rates::rate_N15_p_to_n_O15_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // O15_n_to_He4_C12_reaclib
            reaclib_rates::rate_O15_n_to_He4_C12_reaclib(rate_eval, tfactors);
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
            // O17_n_to_He4_C14_reaclib
            reaclib_rates::rate_O17_n_to_He4_C14_reaclib(rate_eval, tfactors);
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
            // O17_He4_to_n_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O17);
            }
#endif
            reaclib_rates::rate_O17_He4_to_n_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // O18_He4_to_n_Ne21_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O18);
            }
#endif
            reaclib_rates::rate_O18_He4_to_n_Ne21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // F17_n_to_p_O17_reaclib
            reaclib_rates::rate_F17_n_to_p_O17_reaclib(rate_eval, tfactors);
        }

        {
            // F18_n_to_p_O18_reaclib
            reaclib_rates::rate_F18_n_to_p_O18_reaclib(rate_eval, tfactors);
        }

        {
            // F18_n_to_He4_N15_reaclib
            reaclib_rates::rate_F18_n_to_He4_N15_reaclib(rate_eval, tfactors);
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
            // F19_p_to_n_Ne19_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F19);
            }
#endif
            reaclib_rates::rate_F19_p_to_n_Ne19_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ne18_n_to_p_F18_reaclib
            reaclib_rates::rate_Ne18_n_to_p_F18_reaclib(rate_eval, tfactors);
        }

        {
            // Ne18_n_to_He4_O15_reaclib
            reaclib_rates::rate_Ne18_n_to_He4_O15_reaclib(rate_eval, tfactors);
        }

        {
            // Ne18_He4_to_p_Na21_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne18);
            }
#endif
            reaclib_rates::rate_Ne18_He4_to_p_Na21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne19_n_to_He4_O16_reaclib
            reaclib_rates::rate_Ne19_n_to_He4_O16_reaclib(rate_eval, tfactors);
        }

        {
            // Ne19_He4_to_p_Na22_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne19);
            }
#endif
            reaclib_rates::rate_Ne19_He4_to_p_Na22_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_p_to_He4_F17_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne20);
            }
#endif
            reaclib_rates::rate_Ne20_p_to_He4_F17_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ne21_He4_to_n_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne21);
            }
#endif
            reaclib_rates::rate_Ne21_He4_to_n_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne22_He4_to_n_Mg25_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne22);
            }
#endif
            reaclib_rates::rate_Ne22_He4_to_n_Mg25_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na21_n_to_p_Ne21_reaclib
            reaclib_rates::rate_Na21_n_to_p_Ne21_reaclib(rate_eval, tfactors);
        }

        {
            // Na21_n_to_He4_F18_reaclib
            reaclib_rates::rate_Na21_n_to_He4_F18_reaclib(rate_eval, tfactors);
        }

        {
            // Na22_n_to_p_Ne22_reaclib
            reaclib_rates::rate_Na22_n_to_p_Ne22_reaclib(rate_eval, tfactors);
        }

        {
            // Na22_n_to_He4_F19_reaclib
            reaclib_rates::rate_Na22_n_to_He4_F19_reaclib(rate_eval, tfactors);
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
            // Na23_p_to_n_Mg23_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            reaclib_rates::rate_Na23_p_to_n_Mg23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Na23_He4_to_p_Mg26_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na23);
            }
#endif
            reaclib_rates::rate_Na23_He4_to_p_Mg26_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Mg23_He4_to_p_Al26_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg23);
            }
#endif
            reaclib_rates::rate_Mg23_He4_to_p_Al26_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_p_to_He4_Na21_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg24);
            }
#endif
            reaclib_rates::rate_Mg24_p_to_He4_Na21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg25_He4_to_n_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg25);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg25);
            }
#endif
            reaclib_rates::rate_Mg25_He4_to_n_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg26_He4_to_n_Si29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg26);
            }
#endif
            reaclib_rates::rate_Mg26_He4_to_n_Si29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al25_n_to_p_Mg25_reaclib
            reaclib_rates::rate_Al25_n_to_p_Mg25_reaclib(rate_eval, tfactors);
        }

        {
            // Al25_n_to_He4_Na22_reaclib
            reaclib_rates::rate_Al25_n_to_He4_Na22_reaclib(rate_eval, tfactors);
        }

        {
            // Al25_He4_to_p_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al25);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al25);
            }
#endif
            reaclib_rates::rate_Al25_He4_to_p_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al26_n_to_p_Mg26_reaclib
            reaclib_rates::rate_Al26_n_to_p_Mg26_reaclib(rate_eval, tfactors);
        }

        {
            // Al26_n_to_He4_Na23_reaclib
            reaclib_rates::rate_Al26_n_to_He4_Na23_reaclib(rate_eval, tfactors);
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
            // Al27_He4_to_n_P30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al27);
            }
#endif
            reaclib_rates::rate_Al27_He4_to_n_P30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // P29_n_to_He4_Al26_reaclib
            reaclib_rates::rate_P29_n_to_He4_Al26_reaclib(rate_eval, tfactors);
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
            // P30_n_to_p_Si30_reaclib
            reaclib_rates::rate_P30_n_to_p_Si30_reaclib(rate_eval, tfactors);
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
            // P32_n_to_p_Si32_reaclib
            reaclib_rates::rate_P32_n_to_p_Si32_reaclib(rate_eval, tfactors);
        }

        {
            // P32_p_to_n_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P32);
            }
#endif
            reaclib_rates::rate_P32_p_to_n_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P32_p_to_He4_Si29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P32);
            }
#endif
            reaclib_rates::rate_P32_p_to_He4_Si29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // S33_n_to_p_P33_reaclib
            reaclib_rates::rate_S33_n_to_p_P33_reaclib(rate_eval, tfactors);
        }

        {
            // S33_n_to_He4_Si30_reaclib
            reaclib_rates::rate_S33_n_to_He4_Si30_reaclib(rate_eval, tfactors);
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
            // S35_n_to_He4_Si32_reaclib
            reaclib_rates::rate_S35_n_to_He4_Si32_reaclib(rate_eval, tfactors);
        }

        {
            // S35_p_to_He4_P32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S35);
            }
#endif
            reaclib_rates::rate_S35_p_to_He4_P32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // S36_p_to_He4_P33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S36);
            }
#endif
            reaclib_rates::rate_S36_p_to_He4_P33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl33_n_to_p_S33_reaclib
            reaclib_rates::rate_Cl33_n_to_p_S33_reaclib(rate_eval, tfactors);
        }

        {
            // Cl33_n_to_He4_P30_reaclib
            reaclib_rates::rate_Cl33_n_to_He4_P30_reaclib(rate_eval, tfactors);
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
            // Cl34_n_to_p_S34_reaclib
            reaclib_rates::rate_Cl34_n_to_p_S34_reaclib(rate_eval, tfactors);
        }

        {
            // Cl34_n_to_He4_P31_reaclib
            reaclib_rates::rate_Cl34_n_to_He4_P31_reaclib(rate_eval, tfactors);
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
            // Cl35_n_to_p_S35_reaclib
            reaclib_rates::rate_Cl35_n_to_p_S35_reaclib(rate_eval, tfactors);
        }

        {
            // Cl35_n_to_He4_P32_reaclib
            reaclib_rates::rate_Cl35_n_to_He4_P32_reaclib(rate_eval, tfactors);
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
            // Cl36_n_to_p_S36_reaclib
            reaclib_rates::rate_Cl36_n_to_p_S36_reaclib(rate_eval, tfactors);
        }

        {
            // Cl36_n_to_He4_P33_reaclib
            reaclib_rates::rate_Cl36_n_to_He4_P33_reaclib(rate_eval, tfactors);
        }

        {
            // Cl36_p_to_He4_S33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl36);
            }
#endif
            reaclib_rates::rate_Cl36_p_to_He4_S33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar36_n_to_p_Cl36_reaclib
            reaclib_rates::rate_Ar36_n_to_p_Cl36_reaclib(rate_eval, tfactors);
        }

        {
            // Ar36_n_to_He4_S33_reaclib
            reaclib_rates::rate_Ar36_n_to_He4_S33_reaclib(rate_eval, tfactors);
        }

        {
            // Ar37_n_to_p_Cl37_reaclib
            reaclib_rates::rate_Ar37_n_to_p_Cl37_reaclib(rate_eval, tfactors);
        }

        {
            // Ar39_n_to_He4_S36_reaclib
            reaclib_rates::rate_Ar39_n_to_He4_S36_reaclib(rate_eval, tfactors);
        }

        {
            // Ar39_p_to_He4_Cl36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar39);
            }
#endif
            reaclib_rates::rate_Ar39_p_to_He4_Cl36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar40_p_to_He4_Cl37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar40);
            }
#endif
            reaclib_rates::rate_Ar40_p_to_He4_Cl37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K37_n_to_p_Ar37_reaclib
            reaclib_rates::rate_K37_n_to_p_Ar37_reaclib(rate_eval, tfactors);
        }

        {
            // K37_n_to_He4_Cl34_reaclib
            reaclib_rates::rate_K37_n_to_He4_Cl34_reaclib(rate_eval, tfactors);
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
            // K38_n_to_p_Ar38_reaclib
            reaclib_rates::rate_K38_n_to_p_Ar38_reaclib(rate_eval, tfactors);
        }

        {
            // K38_n_to_He4_Cl35_reaclib
            reaclib_rates::rate_K38_n_to_He4_Cl35_reaclib(rate_eval, tfactors);
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
            // K39_n_to_p_Ar39_reaclib
            reaclib_rates::rate_K39_n_to_p_Ar39_reaclib(rate_eval, tfactors);
        }

        {
            // K39_n_to_He4_Cl36_reaclib
            reaclib_rates::rate_K39_n_to_He4_Cl36_reaclib(rate_eval, tfactors);
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
            // K40_n_to_p_Ar40_reaclib
            reaclib_rates::rate_K40_n_to_p_Ar40_reaclib(rate_eval, tfactors);
        }

        {
            // K40_n_to_He4_Cl37_reaclib
            reaclib_rates::rate_K40_n_to_He4_Cl37_reaclib(rate_eval, tfactors);
        }

        {
            // K40_p_to_n_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K40);
            }
#endif
            reaclib_rates::rate_K40_p_to_n_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K40_p_to_He4_Ar37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K40);
            }
#endif
            reaclib_rates::rate_K40_p_to_He4_Ar37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K40_He4_to_p_Ca43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K40);
            }
#endif
            reaclib_rates::rate_K40_He4_to_p_Ca43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K41_p_to_n_Ca41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K41);
            }
#endif
            reaclib_rates::rate_K41_p_to_n_Ca41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K41_p_to_He4_Ar38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K41);
            }
#endif
            reaclib_rates::rate_K41_p_to_He4_Ar38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K41_He4_to_n_Sc44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K41);
            }
#endif
            reaclib_rates::rate_K41_He4_to_n_Sc44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K41_He4_to_p_Ca44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K41);
            }
#endif
            reaclib_rates::rate_K41_He4_to_p_Ca44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_n_to_He4_Ar37_reaclib
            reaclib_rates::rate_Ca40_n_to_He4_Ar37_reaclib(rate_eval, tfactors);
        }

        {
            // Ca41_n_to_He4_Ar38_reaclib
            reaclib_rates::rate_Ca41_n_to_He4_Ar38_reaclib(rate_eval, tfactors);
        }

        {
            // Ca42_n_to_He4_Ar39_reaclib
            reaclib_rates::rate_Ca42_n_to_He4_Ar39_reaclib(rate_eval, tfactors);
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
            // Ca43_n_to_He4_Ar40_reaclib
            reaclib_rates::rate_Ca43_n_to_He4_Ar40_reaclib(rate_eval, tfactors);
        }

        {
            // Ca43_He4_to_n_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca43);
            }
#endif
            reaclib_rates::rate_Ca43_He4_to_n_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca45_He4_to_n_Ti48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca45);
            }
#endif
            reaclib_rates::rate_Ca45_He4_to_n_Ti48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca47_p_to_n_Sc47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca47);
            }
#endif
            reaclib_rates::rate_Ca47_p_to_n_Sc47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca47_He4_to_n_Ti50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca47);
            }
#endif
            reaclib_rates::rate_Ca47_He4_to_n_Ti50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_n_to_p_Ca43_reaclib
            reaclib_rates::rate_Sc43_n_to_p_Ca43_reaclib(rate_eval, tfactors);
        }

        {
            // Sc43_n_to_He4_K40_reaclib
            reaclib_rates::rate_Sc43_n_to_He4_K40_reaclib(rate_eval, tfactors);
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
            // Sc44_n_to_p_Ca44_reaclib
            reaclib_rates::rate_Sc44_n_to_p_Ca44_reaclib(rate_eval, tfactors);
        }

        {
            // Sc44_p_to_He4_Ca41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc44);
            }
#endif
            reaclib_rates::rate_Sc44_p_to_He4_Ca41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc44_He4_to_p_Ti47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc44);
            }
#endif
            reaclib_rates::rate_Sc44_He4_to_p_Ti47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc45_n_to_p_Ca45_reaclib
            reaclib_rates::rate_Sc45_n_to_p_Ca45_reaclib(rate_eval, tfactors);
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
            // Sc46_n_to_p_Ca46_reaclib
            reaclib_rates::rate_Sc46_n_to_p_Ca46_reaclib(rate_eval, tfactors);
        }

        {
            // Sc46_p_to_n_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc46);
            }
#endif
            reaclib_rates::rate_Sc46_p_to_n_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc46_p_to_He4_Ca43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc46);
            }
#endif
            reaclib_rates::rate_Sc46_p_to_He4_Ca43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc46_He4_to_n_V49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc46);
            }
#endif
            reaclib_rates::rate_Sc46_He4_to_n_V49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc46_He4_to_p_Ti49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc46);
            }
#endif
            reaclib_rates::rate_Sc46_He4_to_p_Ti49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc47_p_to_He4_Ca44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc47);
            }
#endif
            reaclib_rates::rate_Sc47_p_to_He4_Ca44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc47_He4_to_p_Ti50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc47);
            }
#endif
            reaclib_rates::rate_Sc47_He4_to_p_Ti50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc48_n_to_p_Ca48_reaclib
            reaclib_rates::rate_Sc48_n_to_p_Ca48_reaclib(rate_eval, tfactors);
        }

        {
            // Sc48_p_to_n_Ti48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc48);
            }
#endif
            reaclib_rates::rate_Sc48_p_to_n_Ti48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc48_p_to_He4_Ca45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc48);
            }
#endif
            reaclib_rates::rate_Sc48_p_to_He4_Ca45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc48_He4_to_n_V51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc48);
            }
#endif
            reaclib_rates::rate_Sc48_He4_to_n_V51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc48_He4_to_p_Ti51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc48);
            }
#endif
            reaclib_rates::rate_Sc48_He4_to_p_Ti51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc49_p_to_n_Ti49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc49);
            }
#endif
            reaclib_rates::rate_Sc49_p_to_n_Ti49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc49_p_to_He4_Ca46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc49);
            }
#endif
            reaclib_rates::rate_Sc49_p_to_He4_Ca46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_n_to_p_Sc44_reaclib
            reaclib_rates::rate_Ti44_n_to_p_Sc44_reaclib(rate_eval, tfactors);
        }

        {
            // Ti44_n_to_He4_Ca41_reaclib
            reaclib_rates::rate_Ti44_n_to_He4_Ca41_reaclib(rate_eval, tfactors);
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
            // Ti45_n_to_p_Sc45_reaclib
            reaclib_rates::rate_Ti45_n_to_p_Sc45_reaclib(rate_eval, tfactors);
        }

        {
            // Ti45_n_to_He4_Ca42_reaclib
            reaclib_rates::rate_Ti45_n_to_He4_Ca42_reaclib(rate_eval, tfactors);
        }

        {
            // Ti45_He4_to_p_V48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti45);
            }
#endif
            reaclib_rates::rate_Ti45_He4_to_p_V48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti47_n_to_p_Sc47_reaclib
            reaclib_rates::rate_Ti47_n_to_p_Sc47_reaclib(rate_eval, tfactors);
        }

        {
            // Ti47_n_to_He4_Ca44_reaclib
            reaclib_rates::rate_Ti47_n_to_He4_Ca44_reaclib(rate_eval, tfactors);
        }

        {
            // Ti49_n_to_He4_Ca46_reaclib
            reaclib_rates::rate_Ti49_n_to_He4_Ca46_reaclib(rate_eval, tfactors);
        }

        {
            // Ti49_He4_to_n_Cr52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti49);
            }
#endif
            reaclib_rates::rate_Ti49_He4_to_n_Cr52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti51_n_to_He4_Ca48_reaclib
            reaclib_rates::rate_Ti51_n_to_He4_Ca48_reaclib(rate_eval, tfactors);
        }

        {
            // Ti51_p_to_n_V51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti51);
            }
#endif
            reaclib_rates::rate_Ti51_p_to_n_V51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti51_He4_to_n_Cr54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti51);
            }
#endif
            reaclib_rates::rate_Ti51_He4_to_n_Cr54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V46_n_to_p_Ti46_reaclib
            reaclib_rates::rate_V46_n_to_p_Ti46_reaclib(rate_eval, tfactors);
        }

        {
            // V46_n_to_He4_Sc43_reaclib
            reaclib_rates::rate_V46_n_to_He4_Sc43_reaclib(rate_eval, tfactors);
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
            // V47_n_to_p_Ti47_reaclib
            reaclib_rates::rate_V47_n_to_p_Ti47_reaclib(rate_eval, tfactors);
        }

        {
            // V47_n_to_He4_Sc44_reaclib
            reaclib_rates::rate_V47_n_to_He4_Sc44_reaclib(rate_eval, tfactors);
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
            // V48_n_to_p_Ti48_reaclib
            reaclib_rates::rate_V48_n_to_p_Ti48_reaclib(rate_eval, tfactors);
        }

        {
            // V48_n_to_He4_Sc45_reaclib
            reaclib_rates::rate_V48_n_to_He4_Sc45_reaclib(rate_eval, tfactors);
        }

        {
            // V48_He4_to_p_Cr51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V48);
            }
#endif
            reaclib_rates::rate_V48_He4_to_p_Cr51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_n_to_p_Ti49_reaclib
            reaclib_rates::rate_V49_n_to_p_Ti49_reaclib(rate_eval, tfactors);
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
            // V50_n_to_p_Ti50_reaclib
            reaclib_rates::rate_V50_n_to_p_Ti50_reaclib(rate_eval, tfactors);
        }

        {
            // V50_n_to_He4_Sc47_reaclib
            reaclib_rates::rate_V50_n_to_He4_Sc47_reaclib(rate_eval, tfactors);
        }

        {
            // V50_p_to_n_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V50);
            }
#endif
            reaclib_rates::rate_V50_p_to_n_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V50_p_to_He4_Ti47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V50);
            }
#endif
            reaclib_rates::rate_V50_p_to_He4_Ti47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V50_He4_to_p_Cr53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V50);
            }
#endif
            reaclib_rates::rate_V50_He4_to_p_Cr53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // V52_n_to_He4_Sc49_reaclib
            reaclib_rates::rate_V52_n_to_He4_Sc49_reaclib(rate_eval, tfactors);
        }

        {
            // V52_p_to_n_Cr52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V52);
            }
#endif
            reaclib_rates::rate_V52_p_to_n_Cr52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V52_p_to_He4_Ti49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V52);
            }
#endif
            reaclib_rates::rate_V52_p_to_He4_Ti49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V52_He4_to_n_Mn55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V52);
            }
#endif
            reaclib_rates::rate_V52_He4_to_n_Mn55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr48_n_to_p_V48_reaclib
            reaclib_rates::rate_Cr48_n_to_p_V48_reaclib(rate_eval, tfactors);
        }

        {
            // Cr48_n_to_He4_Ti45_reaclib
            reaclib_rates::rate_Cr48_n_to_He4_Ti45_reaclib(rate_eval, tfactors);
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
            // Cr49_n_to_p_V49_reaclib
            reaclib_rates::rate_Cr49_n_to_p_V49_reaclib(rate_eval, tfactors);
        }

        {
            // Cr49_n_to_He4_Ti46_reaclib
            reaclib_rates::rate_Cr49_n_to_He4_Ti46_reaclib(rate_eval, tfactors);
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
            // Cr50_n_to_He4_Ti47_reaclib
            reaclib_rates::rate_Cr50_n_to_He4_Ti47_reaclib(rate_eval, tfactors);
        }

        {
            // Cr51_n_to_p_V51_reaclib
            reaclib_rates::rate_Cr51_n_to_p_V51_reaclib(rate_eval, tfactors);
        }

        {
            // Cr51_n_to_He4_Ti48_reaclib
            reaclib_rates::rate_Cr51_n_to_He4_Ti48_reaclib(rate_eval, tfactors);
        }

        {
            // Cr53_n_to_He4_Ti50_reaclib
            reaclib_rates::rate_Cr53_n_to_He4_Ti50_reaclib(rate_eval, tfactors);
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
            // Mn50_n_to_p_Cr50_reaclib
            reaclib_rates::rate_Mn50_n_to_p_Cr50_reaclib(rate_eval, tfactors);
        }

        {
            // Mn50_n_to_He4_V47_reaclib
            reaclib_rates::rate_Mn50_n_to_He4_V47_reaclib(rate_eval, tfactors);
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
            // Mn51_n_to_p_Cr51_reaclib
            reaclib_rates::rate_Mn51_n_to_p_Cr51_reaclib(rate_eval, tfactors);
        }

        {
            // Mn51_n_to_He4_V48_reaclib
            reaclib_rates::rate_Mn51_n_to_He4_V48_reaclib(rate_eval, tfactors);
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
            // Mn52_n_to_p_Cr52_reaclib
            reaclib_rates::rate_Mn52_n_to_p_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // Mn52_n_to_He4_V49_reaclib
            reaclib_rates::rate_Mn52_n_to_He4_V49_reaclib(rate_eval, tfactors);
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
            // Mn53_n_to_p_Cr53_reaclib
            reaclib_rates::rate_Mn53_n_to_p_Cr53_reaclib(rate_eval, tfactors);
        }

        {
            // Mn53_n_to_He4_V50_reaclib
            reaclib_rates::rate_Mn53_n_to_He4_V50_reaclib(rate_eval, tfactors);
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
            // Mn54_n_to_p_Cr54_reaclib
            reaclib_rates::rate_Mn54_n_to_p_Cr54_reaclib(rate_eval, tfactors);
        }

        {
            // Mn54_n_to_He4_V51_reaclib
            reaclib_rates::rate_Mn54_n_to_He4_V51_reaclib(rate_eval, tfactors);
        }

        {
            // Mn54_p_to_He4_Cr51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn54);
            }
#endif
            reaclib_rates::rate_Mn54_p_to_He4_Cr51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Fe52_n_to_p_Mn52_reaclib
            reaclib_rates::rate_Fe52_n_to_p_Mn52_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_n_to_He4_Cr49_reaclib
            reaclib_rates::rate_Fe52_n_to_He4_Cr49_reaclib(rate_eval, tfactors);
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
            // Fe53_n_to_p_Mn53_reaclib
            reaclib_rates::rate_Fe53_n_to_p_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_n_to_He4_Cr50_reaclib
            reaclib_rates::rate_Fe53_n_to_He4_Cr50_reaclib(rate_eval, tfactors);
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
            // Fe54_n_to_p_Mn54_reaclib
            reaclib_rates::rate_Fe54_n_to_p_Mn54_reaclib(rate_eval, tfactors);
        }

        {
            // Fe54_n_to_He4_Cr51_reaclib
            reaclib_rates::rate_Fe54_n_to_He4_Cr51_reaclib(rate_eval, tfactors);
        }

        {
            // Fe55_n_to_p_Mn55_reaclib
            reaclib_rates::rate_Fe55_n_to_p_Mn55_reaclib(rate_eval, tfactors);
        }

        {
            // Fe55_n_to_He4_Cr52_reaclib
            reaclib_rates::rate_Fe55_n_to_He4_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // Fe56_n_to_He4_Cr53_reaclib
            reaclib_rates::rate_Fe56_n_to_He4_Cr53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe57_n_to_He4_Cr54_reaclib
            reaclib_rates::rate_Fe57_n_to_He4_Cr54_reaclib(rate_eval, tfactors);
        }

        {
            // Fe57_p_to_He4_Mn54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe57);
            }
#endif
            reaclib_rates::rate_Fe57_p_to_He4_Mn54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Co53_n_to_p_Fe53_reaclib
            reaclib_rates::rate_Co53_n_to_p_Fe53_reaclib(rate_eval, tfactors);
        }

        {
            // Co53_n_to_He4_Mn50_reaclib
            reaclib_rates::rate_Co53_n_to_He4_Mn50_reaclib(rate_eval, tfactors);
        }

        {
            // Co53_He4_to_p_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co53);
            }
#endif
            reaclib_rates::rate_Co53_He4_to_p_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co54_n_to_p_Fe54_reaclib
            reaclib_rates::rate_Co54_n_to_p_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Co54_n_to_He4_Mn51_reaclib
            reaclib_rates::rate_Co54_n_to_He4_Mn51_reaclib(rate_eval, tfactors);
        }

        {
            // Co54_He4_to_p_Ni57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co54);
            }
#endif
            reaclib_rates::rate_Co54_He4_to_p_Ni57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Co56_n_to_He4_Mn53_reaclib
            reaclib_rates::rate_Co56_n_to_He4_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Co56_He4_to_p_Ni59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co56);
            }
#endif
            reaclib_rates::rate_Co56_He4_to_p_Ni59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co57_n_to_p_Fe57_reaclib
            reaclib_rates::rate_Co57_n_to_p_Fe57_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_n_to_He4_Mn54_reaclib
            reaclib_rates::rate_Co57_n_to_He4_Mn54_reaclib(rate_eval, tfactors);
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
            // Co58_n_to_p_Fe58_reaclib
            reaclib_rates::rate_Co58_n_to_p_Fe58_reaclib(rate_eval, tfactors);
        }

        {
            // Co58_n_to_He4_Mn55_reaclib
            reaclib_rates::rate_Co58_n_to_He4_Mn55_reaclib(rate_eval, tfactors);
        }

        {
            // Co58_p_to_He4_Fe55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co58);
            }
#endif
            reaclib_rates::rate_Co58_p_to_He4_Fe55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Ni58_n_to_p_Co58_reaclib
            reaclib_rates::rate_Ni58_n_to_p_Co58_reaclib(rate_eval, tfactors);
        }

        {
            // Ni58_n_to_He4_Fe55_reaclib
            reaclib_rates::rate_Ni58_n_to_He4_Fe55_reaclib(rate_eval, tfactors);
        }

        {
            // Ni59_n_to_p_Co59_reaclib
            reaclib_rates::rate_Ni59_n_to_p_Co59_reaclib(rate_eval, tfactors);
        }

        {
            // Ni59_n_to_He4_Fe56_reaclib
            reaclib_rates::rate_Ni59_n_to_He4_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Ni60_n_to_He4_Fe57_reaclib
            reaclib_rates::rate_Ni60_n_to_He4_Fe57_reaclib(rate_eval, tfactors);
        }

        {
            // Ni61_n_to_He4_Fe58_reaclib
            reaclib_rates::rate_Ni61_n_to_He4_Fe58_reaclib(rate_eval, tfactors);
        }

        {
            // Ni61_p_to_He4_Co58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni61);
            }
#endif
            reaclib_rates::rate_Ni61_p_to_He4_Co58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni62_p_to_He4_Co59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni62);
            }
#endif
            reaclib_rates::rate_Ni62_p_to_He4_Co59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu57_n_to_p_Ni57_reaclib
            reaclib_rates::rate_Cu57_n_to_p_Ni57_reaclib(rate_eval, tfactors);
        }

        {
            // Cu57_n_to_He4_Co54_reaclib
            reaclib_rates::rate_Cu57_n_to_He4_Co54_reaclib(rate_eval, tfactors);
        }

        {
            // Cu57_He4_to_p_Zn60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu57);
            }
#endif
            reaclib_rates::rate_Cu57_He4_to_p_Zn60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu58_n_to_p_Ni58_reaclib
            reaclib_rates::rate_Cu58_n_to_p_Ni58_reaclib(rate_eval, tfactors);
        }

        {
            // Cu58_n_to_He4_Co55_reaclib
            reaclib_rates::rate_Cu58_n_to_He4_Co55_reaclib(rate_eval, tfactors);
        }

        {
            // Cu59_n_to_p_Ni59_reaclib
            reaclib_rates::rate_Cu59_n_to_p_Ni59_reaclib(rate_eval, tfactors);
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
            // Cu60_n_to_p_Ni60_reaclib
            reaclib_rates::rate_Cu60_n_to_p_Ni60_reaclib(rate_eval, tfactors);
        }

        {
            // Cu60_n_to_He4_Co57_reaclib
            reaclib_rates::rate_Cu60_n_to_He4_Co57_reaclib(rate_eval, tfactors);
        }

        {
            // Cu60_p_to_He4_Ni57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu60);
            }
#endif
            reaclib_rates::rate_Cu60_p_to_He4_Ni57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu61_n_to_p_Ni61_reaclib
            reaclib_rates::rate_Cu61_n_to_p_Ni61_reaclib(rate_eval, tfactors);
        }

        {
            // Cu61_n_to_He4_Co58_reaclib
            reaclib_rates::rate_Cu61_n_to_He4_Co58_reaclib(rate_eval, tfactors);
        }

        {
            // Cu61_p_to_He4_Ni58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu61);
            }
#endif
            reaclib_rates::rate_Cu61_p_to_He4_Ni58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu62_n_to_p_Ni62_reaclib
            reaclib_rates::rate_Cu62_n_to_p_Ni62_reaclib(rate_eval, tfactors);
        }

        {
            // Cu62_n_to_He4_Co59_reaclib
            reaclib_rates::rate_Cu62_n_to_He4_Co59_reaclib(rate_eval, tfactors);
        }

        {
            // Cu62_p_to_He4_Ni59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu62);
            }
#endif
            reaclib_rates::rate_Cu62_p_to_He4_Ni59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu63_n_to_p_Ni63_reaclib
            reaclib_rates::rate_Cu63_n_to_p_Ni63_reaclib(rate_eval, tfactors);
        }

        {
            // Cu63_p_to_He4_Ni60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu63);
            }
#endif
            reaclib_rates::rate_Cu63_p_to_He4_Ni60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu64_n_to_p_Ni64_reaclib
            reaclib_rates::rate_Cu64_n_to_p_Ni64_reaclib(rate_eval, tfactors);
        }

        {
            // Cu64_p_to_He4_Ni61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu64);
            }
#endif
            reaclib_rates::rate_Cu64_p_to_He4_Ni61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cu65_p_to_He4_Ni62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu65);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu65);
            }
#endif
            reaclib_rates::rate_Cu65_p_to_He4_Ni62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn59_n_to_p_Cu59_reaclib
            reaclib_rates::rate_Zn59_n_to_p_Cu59_reaclib(rate_eval, tfactors);
        }

        {
            // Zn59_n_to_He4_Ni56_reaclib
            reaclib_rates::rate_Zn59_n_to_He4_Ni56_reaclib(rate_eval, tfactors);
        }

        {
            // Zn60_n_to_p_Cu60_reaclib
            reaclib_rates::rate_Zn60_n_to_p_Cu60_reaclib(rate_eval, tfactors);
        }

        {
            // Zn60_n_to_He4_Ni57_reaclib
            reaclib_rates::rate_Zn60_n_to_He4_Ni57_reaclib(rate_eval, tfactors);
        }

        {
            // Zn61_n_to_p_Cu61_reaclib
            reaclib_rates::rate_Zn61_n_to_p_Cu61_reaclib(rate_eval, tfactors);
        }

        {
            // Zn61_n_to_He4_Ni58_reaclib
            reaclib_rates::rate_Zn61_n_to_He4_Ni58_reaclib(rate_eval, tfactors);
        }

        {
            // Zn61_p_to_He4_Cu58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn61);
            }
#endif
            reaclib_rates::rate_Zn61_p_to_He4_Cu58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn62_n_to_p_Cu62_reaclib
            reaclib_rates::rate_Zn62_n_to_p_Cu62_reaclib(rate_eval, tfactors);
        }

        {
            // Zn62_n_to_He4_Ni59_reaclib
            reaclib_rates::rate_Zn62_n_to_He4_Ni59_reaclib(rate_eval, tfactors);
        }

        {
            // Zn62_p_to_He4_Cu59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn62);
            }
#endif
            reaclib_rates::rate_Zn62_p_to_He4_Cu59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn63_n_to_p_Cu63_reaclib
            reaclib_rates::rate_Zn63_n_to_p_Cu63_reaclib(rate_eval, tfactors);
        }

        {
            // Zn63_n_to_He4_Ni60_reaclib
            reaclib_rates::rate_Zn63_n_to_He4_Ni60_reaclib(rate_eval, tfactors);
        }

        {
            // Zn63_p_to_He4_Cu60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn63);
            }
#endif
            reaclib_rates::rate_Zn63_p_to_He4_Cu60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn64_n_to_p_Cu64_reaclib
            reaclib_rates::rate_Zn64_n_to_p_Cu64_reaclib(rate_eval, tfactors);
        }

        {
            // Zn64_n_to_He4_Ni61_reaclib
            reaclib_rates::rate_Zn64_n_to_He4_Ni61_reaclib(rate_eval, tfactors);
        }

        {
            // Zn64_p_to_He4_Cu61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn64);
            }
#endif
            reaclib_rates::rate_Zn64_p_to_He4_Cu61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn65_n_to_p_Cu65_reaclib
            reaclib_rates::rate_Zn65_n_to_p_Cu65_reaclib(rate_eval, tfactors);
        }

        {
            // Zn65_n_to_He4_Ni62_reaclib
            reaclib_rates::rate_Zn65_n_to_He4_Ni62_reaclib(rate_eval, tfactors);
        }

        {
            // Zn65_p_to_He4_Cu62_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn65);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn65);
            }
#endif
            reaclib_rates::rate_Zn65_p_to_He4_Cu62_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Zn66_n_to_He4_Ni63_reaclib
            reaclib_rates::rate_Zn66_n_to_He4_Ni63_reaclib(rate_eval, tfactors);
        }

        {
            // Zn66_p_to_He4_Cu63_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn66);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn66);
            }
#endif
            reaclib_rates::rate_Zn66_p_to_He4_Cu63_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ga62_n_to_p_Zn62_reaclib
            reaclib_rates::rate_Ga62_n_to_p_Zn62_reaclib(rate_eval, tfactors);
        }

        {
            // Ga62_n_to_He4_Cu59_reaclib
            reaclib_rates::rate_Ga62_n_to_He4_Cu59_reaclib(rate_eval, tfactors);
        }

        {
            // Ga62_p_to_He4_Zn59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ga62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ga62);
            }
#endif
            reaclib_rates::rate_Ga62_p_to_He4_Zn59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ga63_n_to_p_Zn63_reaclib
            reaclib_rates::rate_Ga63_n_to_p_Zn63_reaclib(rate_eval, tfactors);
        }

        {
            // Ga63_n_to_He4_Cu60_reaclib
            reaclib_rates::rate_Ga63_n_to_He4_Cu60_reaclib(rate_eval, tfactors);
        }

        {
            // Ga63_p_to_He4_Zn60_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ga63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ga63);
            }
#endif
            reaclib_rates::rate_Ga63_p_to_He4_Zn60_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ga64_n_to_p_Zn64_reaclib
            reaclib_rates::rate_Ga64_n_to_p_Zn64_reaclib(rate_eval, tfactors);
        }

        {
            // Ga64_n_to_He4_Cu61_reaclib
            reaclib_rates::rate_Ga64_n_to_He4_Cu61_reaclib(rate_eval, tfactors);
        }

        {
            // Ga64_p_to_He4_Zn61_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ga64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ga64);
            }
#endif
            reaclib_rates::rate_Ga64_p_to_He4_Zn61_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ge63_n_to_p_Ga63_reaclib
            reaclib_rates::rate_Ge63_n_to_p_Ga63_reaclib(rate_eval, tfactors);
        }

        {
            // Ge63_n_to_He4_Zn60_reaclib
            reaclib_rates::rate_Ge63_n_to_He4_Zn60_reaclib(rate_eval, tfactors);
        }

        {
            // Ge64_n_to_p_Ga64_reaclib
            reaclib_rates::rate_Ge64_n_to_p_Ga64_reaclib(rate_eval, tfactors);
        }

        {
            // Ge64_n_to_He4_Zn61_reaclib
            reaclib_rates::rate_Ge64_n_to_He4_Zn61_reaclib(rate_eval, tfactors);
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
            // Li7_d_to_n_He4_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_d_Li7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_d_Li7);
            }
#endif
            reaclib_rates::rate_Li7_d_to_n_He4_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Be9_p_to_d_He4_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Be9);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Be9);
            }
#endif
            reaclib_rates::rate_Be9_p_to_d_He4_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // B8_n_to_p_He4_He4_reaclib
            reaclib_rates::rate_B8_n_to_p_He4_He4_reaclib(rate_eval, tfactors);
        }

        {
            // B11_p_to_3He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_B11);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_B11);
            }
#endif
            reaclib_rates::rate_B11_p_to_3He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Li7_He3_to_n_p_He4_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He3_Li7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He3_Li7);
            }
#endif
            reaclib_rates::rate_Li7_He3_to_n_p_He4_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // Be9_p_to_n_p_He4_He4_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Be9);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Be9);
            }
#endif
            reaclib_rates::rate_Be9_p_to_n_p_He4_He4_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He4_p_n_to_Li6_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_He4);
            }
#endif
            reaclib_rates::rate_He4_p_n_to_Li6_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // He4_He4_n_to_Be9_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4);
            }
#endif
            reaclib_rates::rate_He4_He4_n_to_Be9_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            // p_p_n_to_p_d_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_p);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_p);
            }
#endif
            reaclib_rates::rate_p_p_n_to_p_d_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na21_to_Ne21_electron_capture_weaktab
            tabular_weak_rates::rate_Na21_to_Ne21_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ne21_to_Na21_beta_neg_weaktab
            tabular_weak_rates::rate_Ne21_to_Na21_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Mg23_to_Na23_electron_capture_weaktab
            tabular_weak_rates::rate_Mg23_to_Na23_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Na23_to_Mg23_beta_neg_weaktab
            tabular_weak_rates::rate_Na23_to_Mg23_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Al25_to_Mg25_electron_capture_weaktab
            tabular_weak_rates::rate_Al25_to_Mg25_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mg25_to_Al25_beta_neg_weaktab
            tabular_weak_rates::rate_Mg25_to_Al25_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Al26_to_Mg26_electron_capture_weaktab
            tabular_weak_rates::rate_Al26_to_Mg26_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mg26_to_Al26_beta_neg_weaktab
            tabular_weak_rates::rate_Mg26_to_Al26_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // P30_to_Si30_electron_capture_weaktab
            tabular_weak_rates::rate_P30_to_Si30_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si30_to_P30_beta_neg_weaktab
            tabular_weak_rates::rate_Si30_to_P30_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // P32_to_S32_beta_neg_weaktab
            tabular_weak_rates::rate_P32_to_S32_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P32_to_Si32_electron_capture_weaktab
            tabular_weak_rates::rate_P32_to_Si32_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S32_to_P32_electron_capture_weaktab
            tabular_weak_rates::rate_S32_to_P32_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si32_to_P32_beta_neg_weaktab
            tabular_weak_rates::rate_Si32_to_P32_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Cl34_to_S34_electron_capture_weaktab
            tabular_weak_rates::rate_Cl34_to_S34_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S34_to_Cl34_beta_neg_weaktab
            tabular_weak_rates::rate_S34_to_Cl34_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Ar36_to_Cl36_electron_capture_weaktab
            tabular_weak_rates::rate_Ar36_to_Cl36_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl36_to_Ar36_beta_neg_weaktab
            tabular_weak_rates::rate_Cl36_to_Ar36_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl36_to_S36_electron_capture_weaktab
            tabular_weak_rates::rate_Cl36_to_S36_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S36_to_Cl36_beta_neg_weaktab
            tabular_weak_rates::rate_S36_to_Cl36_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Ar38_to_K38_beta_neg_weaktab
            tabular_weak_rates::rate_Ar38_to_K38_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K38_to_Ar38_electron_capture_weaktab
            tabular_weak_rates::rate_K38_to_Ar38_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar39_to_K39_beta_neg_weaktab
            tabular_weak_rates::rate_Ar39_to_K39_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K39_to_Ar39_electron_capture_weaktab
            tabular_weak_rates::rate_K39_to_Ar39_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar40_to_K40_beta_neg_weaktab
            tabular_weak_rates::rate_Ar40_to_K40_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca40_to_K40_electron_capture_weaktab
            tabular_weak_rates::rate_Ca40_to_K40_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K40_to_Ar40_electron_capture_weaktab
            tabular_weak_rates::rate_K40_to_Ar40_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K40_to_Ca40_beta_neg_weaktab
            tabular_weak_rates::rate_K40_to_Ca40_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca41_to_K41_electron_capture_weaktab
            tabular_weak_rates::rate_Ca41_to_K41_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K41_to_Ca41_beta_neg_weaktab
            tabular_weak_rates::rate_K41_to_Ca41_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca43_to_Sc43_beta_neg_weaktab
            tabular_weak_rates::rate_Ca43_to_Sc43_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc43_to_Ca43_electron_capture_weaktab
            tabular_weak_rates::rate_Sc43_to_Ca43_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca44_to_Sc44_beta_neg_weaktab
            tabular_weak_rates::rate_Ca44_to_Sc44_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc44_to_Ca44_electron_capture_weaktab
            tabular_weak_rates::rate_Sc44_to_Ca44_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc44_to_Ti44_beta_neg_weaktab
            tabular_weak_rates::rate_Sc44_to_Ti44_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti44_to_Sc44_electron_capture_weaktab
            tabular_weak_rates::rate_Ti44_to_Sc44_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co53_to_Fe53_electron_capture_weaktab
            tabular_weak_rates::rate_Co53_to_Fe53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe53_to_Co53_beta_neg_weaktab
            tabular_weak_rates::rate_Fe53_to_Co53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu57_to_Ni57_electron_capture_weaktab
            tabular_weak_rates::rate_Cu57_to_Ni57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Cu57_beta_neg_weaktab
            tabular_weak_rates::rate_Ni57_to_Cu57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca45_to_Sc45_beta_neg_weaktab
            tabular_weak_rates::rate_Ca45_to_Sc45_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc45_to_Ca45_electron_capture_weaktab
            tabular_weak_rates::rate_Sc45_to_Ca45_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc45_to_Ti45_beta_neg_weaktab
            tabular_weak_rates::rate_Sc45_to_Ti45_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti45_to_Sc45_electron_capture_weaktab
            tabular_weak_rates::rate_Ti45_to_Sc45_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca46_to_Sc46_beta_neg_weaktab
            tabular_weak_rates::rate_Ca46_to_Sc46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc46_to_Ca46_electron_capture_weaktab
            tabular_weak_rates::rate_Sc46_to_Ca46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc46_to_Ti46_beta_neg_weaktab
            tabular_weak_rates::rate_Sc46_to_Ti46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti46_to_Sc46_electron_capture_weaktab
            tabular_weak_rates::rate_Ti46_to_Sc46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti46_to_V46_beta_neg_weaktab
            tabular_weak_rates::rate_Ti46_to_V46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V46_to_Ti46_electron_capture_weaktab
            tabular_weak_rates::rate_V46_to_Ti46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca47_to_Sc47_beta_neg_weaktab
            tabular_weak_rates::rate_Ca47_to_Sc47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc47_to_Ca47_electron_capture_weaktab
            tabular_weak_rates::rate_Sc47_to_Ca47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc47_to_Ti47_beta_neg_weaktab
            tabular_weak_rates::rate_Sc47_to_Ti47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti47_to_Sc47_electron_capture_weaktab
            tabular_weak_rates::rate_Ti47_to_Sc47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti47_to_V47_beta_neg_weaktab
            tabular_weak_rates::rate_Ti47_to_V47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V47_to_Ti47_electron_capture_weaktab
            tabular_weak_rates::rate_V47_to_Ti47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca48_to_Sc48_beta_neg_weaktab
            tabular_weak_rates::rate_Ca48_to_Sc48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr48_to_V48_electron_capture_weaktab
            tabular_weak_rates::rate_Cr48_to_V48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc48_to_Ca48_electron_capture_weaktab
            tabular_weak_rates::rate_Sc48_to_Ca48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc48_to_Ti48_beta_neg_weaktab
            tabular_weak_rates::rate_Sc48_to_Ti48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti48_to_Sc48_electron_capture_weaktab
            tabular_weak_rates::rate_Ti48_to_Sc48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti48_to_V48_beta_neg_weaktab
            tabular_weak_rates::rate_Ti48_to_V48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V48_to_Cr48_beta_neg_weaktab
            tabular_weak_rates::rate_V48_to_Cr48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V48_to_Ti48_electron_capture_weaktab
            tabular_weak_rates::rate_V48_to_Ti48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr49_to_V49_electron_capture_weaktab
            tabular_weak_rates::rate_Cr49_to_V49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc49_to_Ti49_beta_neg_weaktab
            tabular_weak_rates::rate_Sc49_to_Ti49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti49_to_Sc49_electron_capture_weaktab
            tabular_weak_rates::rate_Ti49_to_Sc49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti49_to_V49_beta_neg_weaktab
            tabular_weak_rates::rate_Ti49_to_V49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V49_to_Cr49_beta_neg_weaktab
            tabular_weak_rates::rate_V49_to_Cr49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V49_to_Ti49_electron_capture_weaktab
            tabular_weak_rates::rate_V49_to_Ti49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr50_to_Mn50_beta_neg_weaktab
            tabular_weak_rates::rate_Cr50_to_Mn50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr50_to_V50_electron_capture_weaktab
            tabular_weak_rates::rate_Cr50_to_V50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn50_to_Cr50_electron_capture_weaktab
            tabular_weak_rates::rate_Mn50_to_Cr50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti50_to_V50_beta_neg_weaktab
            tabular_weak_rates::rate_Ti50_to_V50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V50_to_Cr50_beta_neg_weaktab
            tabular_weak_rates::rate_V50_to_Cr50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V50_to_Ti50_electron_capture_weaktab
            tabular_weak_rates::rate_V50_to_Ti50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr51_to_Mn51_beta_neg_weaktab
            tabular_weak_rates::rate_Cr51_to_Mn51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr51_to_V51_electron_capture_weaktab
            tabular_weak_rates::rate_Cr51_to_V51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn51_to_Cr51_electron_capture_weaktab
            tabular_weak_rates::rate_Mn51_to_Cr51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti51_to_V51_beta_neg_weaktab
            tabular_weak_rates::rate_Ti51_to_V51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V51_to_Cr51_beta_neg_weaktab
            tabular_weak_rates::rate_V51_to_Cr51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V51_to_Ti51_electron_capture_weaktab
            tabular_weak_rates::rate_V51_to_Ti51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr52_to_Mn52_beta_neg_weaktab
            tabular_weak_rates::rate_Cr52_to_Mn52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr52_to_V52_electron_capture_weaktab
            tabular_weak_rates::rate_Cr52_to_V52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe52_to_Mn52_electron_capture_weaktab
            tabular_weak_rates::rate_Fe52_to_Mn52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Cr52_electron_capture_weaktab
            tabular_weak_rates::rate_Mn52_to_Cr52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Fe52_beta_neg_weaktab
            tabular_weak_rates::rate_Mn52_to_Fe52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V52_to_Cr52_beta_neg_weaktab
            tabular_weak_rates::rate_V52_to_Cr52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr53_to_Mn53_beta_neg_weaktab
            tabular_weak_rates::rate_Cr53_to_Mn53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe53_to_Mn53_electron_capture_weaktab
            tabular_weak_rates::rate_Fe53_to_Mn53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn53_to_Cr53_electron_capture_weaktab
            tabular_weak_rates::rate_Mn53_to_Cr53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn53_to_Fe53_beta_neg_weaktab
            tabular_weak_rates::rate_Mn53_to_Fe53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co54_to_Fe54_electron_capture_weaktab
            tabular_weak_rates::rate_Co54_to_Fe54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr54_to_Mn54_beta_neg_weaktab
            tabular_weak_rates::rate_Cr54_to_Mn54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe54_to_Co54_beta_neg_weaktab
            tabular_weak_rates::rate_Fe54_to_Co54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe54_to_Mn54_electron_capture_weaktab
            tabular_weak_rates::rate_Fe54_to_Mn54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn54_to_Cr54_electron_capture_weaktab
            tabular_weak_rates::rate_Mn54_to_Cr54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn54_to_Fe54_beta_neg_weaktab
            tabular_weak_rates::rate_Mn54_to_Fe54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Fe55_to_Mn55_electron_capture_weaktab
            tabular_weak_rates::rate_Fe55_to_Mn55_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn55_to_Fe55_beta_neg_weaktab
            tabular_weak_rates::rate_Mn55_to_Fe55_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Co57_to_Fe57_electron_capture_weaktab
            tabular_weak_rates::rate_Co57_to_Fe57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co57_to_Ni57_beta_neg_weaktab
            tabular_weak_rates::rate_Co57_to_Ni57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe57_to_Co57_beta_neg_weaktab
            tabular_weak_rates::rate_Fe57_to_Co57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Co57_electron_capture_weaktab
            tabular_weak_rates::rate_Ni57_to_Co57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co58_to_Fe58_electron_capture_weaktab
            tabular_weak_rates::rate_Co58_to_Fe58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co58_to_Ni58_beta_neg_weaktab
            tabular_weak_rates::rate_Co58_to_Ni58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu58_to_Ni58_electron_capture_weaktab
            tabular_weak_rates::rate_Cu58_to_Ni58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe58_to_Co58_beta_neg_weaktab
            tabular_weak_rates::rate_Fe58_to_Co58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni58_to_Co58_electron_capture_weaktab
            tabular_weak_rates::rate_Ni58_to_Co58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni58_to_Cu58_beta_neg_weaktab
            tabular_weak_rates::rate_Ni58_to_Cu58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co59_to_Ni59_beta_neg_weaktab
            tabular_weak_rates::rate_Co59_to_Ni59_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu59_to_Ni59_electron_capture_weaktab
            tabular_weak_rates::rate_Cu59_to_Ni59_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni59_to_Co59_electron_capture_weaktab
            tabular_weak_rates::rate_Ni59_to_Co59_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni59_to_Cu59_beta_neg_weaktab
            tabular_weak_rates::rate_Ni59_to_Cu59_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu60_to_Ni60_electron_capture_weaktab
            tabular_weak_rates::rate_Cu60_to_Ni60_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu60_to_Zn60_beta_neg_weaktab
            tabular_weak_rates::rate_Cu60_to_Zn60_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni60_to_Cu60_beta_neg_weaktab
            tabular_weak_rates::rate_Ni60_to_Cu60_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn60_to_Cu60_electron_capture_weaktab
            tabular_weak_rates::rate_Zn60_to_Cu60_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu61_to_Ni61_electron_capture_weaktab
            tabular_weak_rates::rate_Cu61_to_Ni61_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu61_to_Zn61_beta_neg_weaktab
            tabular_weak_rates::rate_Cu61_to_Zn61_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni61_to_Cu61_beta_neg_weaktab
            tabular_weak_rates::rate_Ni61_to_Cu61_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn61_to_Cu61_electron_capture_weaktab
            tabular_weak_rates::rate_Zn61_to_Cu61_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu62_to_Ni62_electron_capture_weaktab
            tabular_weak_rates::rate_Cu62_to_Ni62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu62_to_Zn62_beta_neg_weaktab
            tabular_weak_rates::rate_Cu62_to_Zn62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga62_to_Zn62_electron_capture_weaktab
            tabular_weak_rates::rate_Ga62_to_Zn62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni62_to_Cu62_beta_neg_weaktab
            tabular_weak_rates::rate_Ni62_to_Cu62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn62_to_Cu62_electron_capture_weaktab
            tabular_weak_rates::rate_Zn62_to_Cu62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn62_to_Ga62_beta_neg_weaktab
            tabular_weak_rates::rate_Zn62_to_Ga62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu63_to_Ni63_electron_capture_weaktab
            tabular_weak_rates::rate_Cu63_to_Ni63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu63_to_Zn63_beta_neg_weaktab
            tabular_weak_rates::rate_Cu63_to_Zn63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga63_to_Zn63_electron_capture_weaktab
            tabular_weak_rates::rate_Ga63_to_Zn63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni63_to_Cu63_beta_neg_weaktab
            tabular_weak_rates::rate_Ni63_to_Cu63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn63_to_Cu63_electron_capture_weaktab
            tabular_weak_rates::rate_Zn63_to_Cu63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn63_to_Ga63_beta_neg_weaktab
            tabular_weak_rates::rate_Zn63_to_Ga63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu64_to_Ni64_electron_capture_weaktab
            tabular_weak_rates::rate_Cu64_to_Ni64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu64_to_Zn64_beta_neg_weaktab
            tabular_weak_rates::rate_Cu64_to_Zn64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga64_to_Ge64_beta_neg_weaktab
            tabular_weak_rates::rate_Ga64_to_Ge64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga64_to_Zn64_electron_capture_weaktab
            tabular_weak_rates::rate_Ga64_to_Zn64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ge64_to_Ga64_electron_capture_weaktab
            tabular_weak_rates::rate_Ge64_to_Ga64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni64_to_Cu64_beta_neg_weaktab
            tabular_weak_rates::rate_Ni64_to_Cu64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn64_to_Cu64_electron_capture_weaktab
            tabular_weak_rates::rate_Zn64_to_Cu64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn64_to_Ga64_beta_neg_weaktab
            tabular_weak_rates::rate_Zn64_to_Ga64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu65_to_Zn65_beta_neg_weaktab
            tabular_weak_rates::rate_Cu65_to_Zn65_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn65_to_Cu65_electron_capture_weaktab
            tabular_weak_rates::rate_Zn65_to_Cu65_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // d_to_n_p_derived
            derived_rates::rate_d_to_n_p_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // He3_to_p_d_derived
            derived_rates::rate_He3_to_p_d_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // He4_to_n_He3_derived
            derived_rates::rate_He4_to_n_He3_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // He4_to_d_d_derived
            derived_rates::rate_He4_to_d_d_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Li6_to_He4_d_derived
            derived_rates::rate_Li6_to_He4_d_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Li7_to_n_Li6_derived
            derived_rates::rate_Li7_to_n_Li6_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Be7_to_p_Li6_derived
            derived_rates::rate_Be7_to_p_Li6_derived(rate_eval, tfactors, pf_cache);
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
            // B10_to_p_Be9_derived
            derived_rates::rate_B10_to_p_Be9_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // B10_to_He4_Li6_derived
            derived_rates::rate_B10_to_He4_Li6_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // B11_to_n_B10_derived
            derived_rates::rate_B11_to_n_B10_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // B11_to_He4_Li7_derived
            derived_rates::rate_B11_to_He4_Li7_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C12_to_p_B11_derived
            derived_rates::rate_C12_to_p_B11_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C13_to_n_C12_derived
            derived_rates::rate_C13_to_n_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C14_to_n_C13_derived
            derived_rates::rate_C14_to_n_C13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N13_to_p_C12_derived
            derived_rates::rate_N13_to_p_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N14_to_n_N13_derived
            derived_rates::rate_N14_to_n_N13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N14_to_p_C13_derived
            derived_rates::rate_N14_to_p_C13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N15_to_n_N14_derived
            derived_rates::rate_N15_to_n_N14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N15_to_p_C14_derived
            derived_rates::rate_N15_to_p_C14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O14_to_p_N13_derived
            derived_rates::rate_O14_to_p_N13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O15_to_n_O14_derived
            derived_rates::rate_O15_to_n_O14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O15_to_p_N14_derived
            derived_rates::rate_O15_to_p_N14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O16_to_n_O15_derived
            derived_rates::rate_O16_to_n_O15_derived(rate_eval, tfactors, pf_cache);
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
            // O17_to_n_O16_derived
            derived_rates::rate_O17_to_n_O16_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O18_to_n_O17_derived
            derived_rates::rate_O18_to_n_O17_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O18_to_He4_C14_derived
            derived_rates::rate_O18_to_He4_C14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // F17_to_p_O16_derived
            derived_rates::rate_F17_to_p_O16_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // F18_to_n_F17_derived
            derived_rates::rate_F18_to_n_F17_derived(rate_eval, tfactors, pf_cache);
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
            // F19_to_n_F18_derived
            derived_rates::rate_F19_to_n_F18_derived(rate_eval, tfactors, pf_cache);
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
            // Ne18_to_p_F17_derived
            derived_rates::rate_Ne18_to_p_F17_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne18_to_He4_O14_derived
            derived_rates::rate_Ne18_to_He4_O14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne19_to_n_Ne18_derived
            derived_rates::rate_Ne19_to_n_Ne18_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne19_to_p_F18_derived
            derived_rates::rate_Ne19_to_p_F18_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne19_to_He4_O15_derived
            derived_rates::rate_Ne19_to_He4_O15_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_to_n_Ne19_derived
            derived_rates::rate_Ne20_to_n_Ne19_derived(rate_eval, tfactors, pf_cache);
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
            // Ne21_to_n_Ne20_derived
            derived_rates::rate_Ne21_to_n_Ne20_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne21_to_He4_O17_derived
            derived_rates::rate_Ne21_to_He4_O17_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne22_to_n_Ne21_derived
            derived_rates::rate_Ne22_to_n_Ne21_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne22_to_He4_O18_derived
            derived_rates::rate_Ne22_to_He4_O18_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na21_to_p_Ne20_derived
            derived_rates::rate_Na21_to_p_Ne20_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na21_to_He4_F17_derived
            derived_rates::rate_Na21_to_He4_F17_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na22_to_n_Na21_derived
            derived_rates::rate_Na22_to_n_Na21_derived(rate_eval, tfactors, pf_cache);
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
            // Na23_to_n_Na22_derived
            derived_rates::rate_Na23_to_n_Na22_derived(rate_eval, tfactors, pf_cache);
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
            // Mg23_to_p_Na22_derived
            derived_rates::rate_Mg23_to_p_Na22_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg23_to_He4_Ne19_derived
            derived_rates::rate_Mg23_to_He4_Ne19_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg24_to_n_Mg23_derived
            derived_rates::rate_Mg24_to_n_Mg23_derived(rate_eval, tfactors, pf_cache);
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
            // Mg25_to_n_Mg24_derived
            derived_rates::rate_Mg25_to_n_Mg24_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg25_to_He4_Ne21_derived
            derived_rates::rate_Mg25_to_He4_Ne21_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg26_to_n_Mg25_derived
            derived_rates::rate_Mg26_to_n_Mg25_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg26_to_He4_Ne22_derived
            derived_rates::rate_Mg26_to_He4_Ne22_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Al25_to_p_Mg24_derived
            derived_rates::rate_Al25_to_p_Mg24_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Al25_to_He4_Na21_derived
            derived_rates::rate_Al25_to_He4_Na21_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Al26_to_n_Al25_derived
            derived_rates::rate_Al26_to_n_Al25_derived(rate_eval, tfactors, pf_cache);
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
            // Al27_to_n_Al26_derived
            derived_rates::rate_Al27_to_n_Al26_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Al27_to_p_Mg26_derived
            derived_rates::rate_Al27_to_p_Mg26_derived(rate_eval, tfactors, pf_cache);
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
            // Si29_to_n_Si28_derived
            derived_rates::rate_Si29_to_n_Si28_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si29_to_He4_Mg25_derived
            derived_rates::rate_Si29_to_He4_Mg25_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si30_to_n_Si29_derived
            derived_rates::rate_Si30_to_n_Si29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si30_to_He4_Mg26_derived
            derived_rates::rate_Si30_to_He4_Mg26_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si31_to_n_Si30_derived
            derived_rates::rate_Si31_to_n_Si30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si32_to_n_Si31_derived
            derived_rates::rate_Si32_to_n_Si31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P29_to_p_Si28_derived
            derived_rates::rate_P29_to_p_Si28_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P29_to_He4_Al25_derived
            derived_rates::rate_P29_to_He4_Al25_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P30_to_n_P29_derived
            derived_rates::rate_P30_to_n_P29_derived(rate_eval, tfactors, pf_cache);
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
            // P31_to_n_P30_derived
            derived_rates::rate_P31_to_n_P30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P31_to_p_Si30_derived
            derived_rates::rate_P31_to_p_Si30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P31_to_He4_Al27_derived
            derived_rates::rate_P31_to_He4_Al27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P32_to_n_P31_derived
            derived_rates::rate_P32_to_n_P31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P32_to_p_Si31_derived
            derived_rates::rate_P32_to_p_Si31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P33_to_n_P32_derived
            derived_rates::rate_P33_to_n_P32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P33_to_p_Si32_derived
            derived_rates::rate_P33_to_p_Si32_derived(rate_eval, tfactors, pf_cache);
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
            // S33_to_n_S32_derived
            derived_rates::rate_S33_to_n_S32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S33_to_p_P32_derived
            derived_rates::rate_S33_to_p_P32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S33_to_He4_Si29_derived
            derived_rates::rate_S33_to_He4_Si29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S34_to_n_S33_derived
            derived_rates::rate_S34_to_n_S33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S34_to_p_P33_derived
            derived_rates::rate_S34_to_p_P33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S34_to_He4_Si30_derived
            derived_rates::rate_S34_to_He4_Si30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S35_to_n_S34_derived
            derived_rates::rate_S35_to_n_S34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S35_to_He4_Si31_derived
            derived_rates::rate_S35_to_He4_Si31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S36_to_n_S35_derived
            derived_rates::rate_S36_to_n_S35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S36_to_He4_Si32_derived
            derived_rates::rate_S36_to_He4_Si32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl33_to_p_S32_derived
            derived_rates::rate_Cl33_to_p_S32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl33_to_He4_P29_derived
            derived_rates::rate_Cl33_to_He4_P29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl34_to_n_Cl33_derived
            derived_rates::rate_Cl34_to_n_Cl33_derived(rate_eval, tfactors, pf_cache);
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
            // Cl35_to_n_Cl34_derived
            derived_rates::rate_Cl35_to_n_Cl34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl35_to_p_S34_derived
            derived_rates::rate_Cl35_to_p_S34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl35_to_He4_P31_derived
            derived_rates::rate_Cl35_to_He4_P31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl36_to_n_Cl35_derived
            derived_rates::rate_Cl36_to_n_Cl35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl36_to_p_S35_derived
            derived_rates::rate_Cl36_to_p_S35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl36_to_He4_P32_derived
            derived_rates::rate_Cl36_to_He4_P32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl37_to_n_Cl36_derived
            derived_rates::rate_Cl37_to_n_Cl36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl37_to_p_S36_derived
            derived_rates::rate_Cl37_to_p_S36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl37_to_He4_P33_derived
            derived_rates::rate_Cl37_to_He4_P33_derived(rate_eval, tfactors, pf_cache);
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
            // Ar37_to_n_Ar36_derived
            derived_rates::rate_Ar37_to_n_Ar36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar37_to_p_Cl36_derived
            derived_rates::rate_Ar37_to_p_Cl36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar37_to_He4_S33_derived
            derived_rates::rate_Ar37_to_He4_S33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar38_to_n_Ar37_derived
            derived_rates::rate_Ar38_to_n_Ar37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar38_to_p_Cl37_derived
            derived_rates::rate_Ar38_to_p_Cl37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar38_to_He4_S34_derived
            derived_rates::rate_Ar38_to_He4_S34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar39_to_n_Ar38_derived
            derived_rates::rate_Ar39_to_n_Ar38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar39_to_He4_S35_derived
            derived_rates::rate_Ar39_to_He4_S35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar40_to_n_Ar39_derived
            derived_rates::rate_Ar40_to_n_Ar39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar40_to_He4_S36_derived
            derived_rates::rate_Ar40_to_He4_S36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K37_to_p_Ar36_derived
            derived_rates::rate_K37_to_p_Ar36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K37_to_He4_Cl33_derived
            derived_rates::rate_K37_to_He4_Cl33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K38_to_n_K37_derived
            derived_rates::rate_K38_to_n_K37_derived(rate_eval, tfactors, pf_cache);
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
            // K39_to_n_K38_derived
            derived_rates::rate_K39_to_n_K38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K39_to_p_Ar38_derived
            derived_rates::rate_K39_to_p_Ar38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K39_to_He4_Cl35_derived
            derived_rates::rate_K39_to_He4_Cl35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K40_to_n_K39_derived
            derived_rates::rate_K40_to_n_K39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K40_to_p_Ar39_derived
            derived_rates::rate_K40_to_p_Ar39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K40_to_He4_Cl36_derived
            derived_rates::rate_K40_to_He4_Cl36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K41_to_n_K40_derived
            derived_rates::rate_K41_to_n_K40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K41_to_p_Ar40_derived
            derived_rates::rate_K41_to_p_Ar40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K41_to_He4_Cl37_derived
            derived_rates::rate_K41_to_He4_Cl37_derived(rate_eval, tfactors, pf_cache);
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
            // Ca41_to_n_Ca40_derived
            derived_rates::rate_Ca41_to_n_Ca40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca41_to_p_K40_derived
            derived_rates::rate_Ca41_to_p_K40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca41_to_He4_Ar37_derived
            derived_rates::rate_Ca41_to_He4_Ar37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca42_to_n_Ca41_derived
            derived_rates::rate_Ca42_to_n_Ca41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca42_to_p_K41_derived
            derived_rates::rate_Ca42_to_p_K41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca42_to_He4_Ar38_derived
            derived_rates::rate_Ca42_to_He4_Ar38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca43_to_n_Ca42_derived
            derived_rates::rate_Ca43_to_n_Ca42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca43_to_He4_Ar39_derived
            derived_rates::rate_Ca43_to_He4_Ar39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca44_to_n_Ca43_derived
            derived_rates::rate_Ca44_to_n_Ca43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca44_to_He4_Ar40_derived
            derived_rates::rate_Ca44_to_He4_Ar40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca45_to_n_Ca44_derived
            derived_rates::rate_Ca45_to_n_Ca44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca46_to_n_Ca45_derived
            derived_rates::rate_Ca46_to_n_Ca45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca47_to_n_Ca46_derived
            derived_rates::rate_Ca47_to_n_Ca46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca48_to_n_Ca47_derived
            derived_rates::rate_Ca48_to_n_Ca47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc43_to_p_Ca42_derived
            derived_rates::rate_Sc43_to_p_Ca42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc43_to_He4_K39_derived
            derived_rates::rate_Sc43_to_He4_K39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc44_to_n_Sc43_derived
            derived_rates::rate_Sc44_to_n_Sc43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc44_to_p_Ca43_derived
            derived_rates::rate_Sc44_to_p_Ca43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc44_to_He4_K40_derived
            derived_rates::rate_Sc44_to_He4_K40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc45_to_n_Sc44_derived
            derived_rates::rate_Sc45_to_n_Sc44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc45_to_p_Ca44_derived
            derived_rates::rate_Sc45_to_p_Ca44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc45_to_He4_K41_derived
            derived_rates::rate_Sc45_to_He4_K41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc46_to_n_Sc45_derived
            derived_rates::rate_Sc46_to_n_Sc45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc46_to_p_Ca45_derived
            derived_rates::rate_Sc46_to_p_Ca45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc47_to_n_Sc46_derived
            derived_rates::rate_Sc47_to_n_Sc46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc47_to_p_Ca46_derived
            derived_rates::rate_Sc47_to_p_Ca46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc48_to_n_Sc47_derived
            derived_rates::rate_Sc48_to_n_Sc47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc48_to_p_Ca47_derived
            derived_rates::rate_Sc48_to_p_Ca47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc49_to_n_Sc48_derived
            derived_rates::rate_Sc49_to_n_Sc48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc49_to_p_Ca48_derived
            derived_rates::rate_Sc49_to_p_Ca48_derived(rate_eval, tfactors, pf_cache);
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
            // Ti45_to_n_Ti44_derived
            derived_rates::rate_Ti45_to_n_Ti44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti45_to_p_Sc44_derived
            derived_rates::rate_Ti45_to_p_Sc44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti45_to_He4_Ca41_derived
            derived_rates::rate_Ti45_to_He4_Ca41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti46_to_n_Ti45_derived
            derived_rates::rate_Ti46_to_n_Ti45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti46_to_p_Sc45_derived
            derived_rates::rate_Ti46_to_p_Sc45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti46_to_He4_Ca42_derived
            derived_rates::rate_Ti46_to_He4_Ca42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti47_to_n_Ti46_derived
            derived_rates::rate_Ti47_to_n_Ti46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti47_to_p_Sc46_derived
            derived_rates::rate_Ti47_to_p_Sc46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti47_to_He4_Ca43_derived
            derived_rates::rate_Ti47_to_He4_Ca43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti48_to_n_Ti47_derived
            derived_rates::rate_Ti48_to_n_Ti47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti48_to_p_Sc47_derived
            derived_rates::rate_Ti48_to_p_Sc47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti48_to_He4_Ca44_derived
            derived_rates::rate_Ti48_to_He4_Ca44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti49_to_n_Ti48_derived
            derived_rates::rate_Ti49_to_n_Ti48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti49_to_p_Sc48_derived
            derived_rates::rate_Ti49_to_p_Sc48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti49_to_He4_Ca45_derived
            derived_rates::rate_Ti49_to_He4_Ca45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti50_to_n_Ti49_derived
            derived_rates::rate_Ti50_to_n_Ti49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti50_to_p_Sc49_derived
            derived_rates::rate_Ti50_to_p_Sc49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti50_to_He4_Ca46_derived
            derived_rates::rate_Ti50_to_He4_Ca46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti51_to_n_Ti50_derived
            derived_rates::rate_Ti51_to_n_Ti50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti51_to_He4_Ca47_derived
            derived_rates::rate_Ti51_to_He4_Ca47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V46_to_p_Ti45_derived
            derived_rates::rate_V46_to_p_Ti45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V47_to_n_V46_derived
            derived_rates::rate_V47_to_n_V46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V47_to_p_Ti46_derived
            derived_rates::rate_V47_to_p_Ti46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V47_to_He4_Sc43_derived
            derived_rates::rate_V47_to_He4_Sc43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V48_to_n_V47_derived
            derived_rates::rate_V48_to_n_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V48_to_p_Ti47_derived
            derived_rates::rate_V48_to_p_Ti47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V48_to_He4_Sc44_derived
            derived_rates::rate_V48_to_He4_Sc44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V49_to_n_V48_derived
            derived_rates::rate_V49_to_n_V48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V49_to_p_Ti48_derived
            derived_rates::rate_V49_to_p_Ti48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V49_to_He4_Sc45_derived
            derived_rates::rate_V49_to_He4_Sc45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V50_to_n_V49_derived
            derived_rates::rate_V50_to_n_V49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V50_to_p_Ti49_derived
            derived_rates::rate_V50_to_p_Ti49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V50_to_He4_Sc46_derived
            derived_rates::rate_V50_to_He4_Sc46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V51_to_n_V50_derived
            derived_rates::rate_V51_to_n_V50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V51_to_p_Ti50_derived
            derived_rates::rate_V51_to_p_Ti50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V51_to_He4_Sc47_derived
            derived_rates::rate_V51_to_He4_Sc47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V52_to_n_V51_derived
            derived_rates::rate_V52_to_n_V51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V52_to_p_Ti51_derived
            derived_rates::rate_V52_to_p_Ti51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V52_to_He4_Sc48_derived
            derived_rates::rate_V52_to_He4_Sc48_derived(rate_eval, tfactors, pf_cache);
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
            // Cr49_to_n_Cr48_derived
            derived_rates::rate_Cr49_to_n_Cr48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr49_to_p_V48_derived
            derived_rates::rate_Cr49_to_p_V48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr49_to_He4_Ti45_derived
            derived_rates::rate_Cr49_to_He4_Ti45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr50_to_n_Cr49_derived
            derived_rates::rate_Cr50_to_n_Cr49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr50_to_p_V49_derived
            derived_rates::rate_Cr50_to_p_V49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr50_to_He4_Ti46_derived
            derived_rates::rate_Cr50_to_He4_Ti46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr51_to_n_Cr50_derived
            derived_rates::rate_Cr51_to_n_Cr50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr51_to_p_V50_derived
            derived_rates::rate_Cr51_to_p_V50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr51_to_He4_Ti47_derived
            derived_rates::rate_Cr51_to_He4_Ti47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr52_to_n_Cr51_derived
            derived_rates::rate_Cr52_to_n_Cr51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr52_to_p_V51_derived
            derived_rates::rate_Cr52_to_p_V51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr52_to_He4_Ti48_derived
            derived_rates::rate_Cr52_to_He4_Ti48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr53_to_n_Cr52_derived
            derived_rates::rate_Cr53_to_n_Cr52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr53_to_p_V52_derived
            derived_rates::rate_Cr53_to_p_V52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr53_to_He4_Ti49_derived
            derived_rates::rate_Cr53_to_He4_Ti49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr54_to_n_Cr53_derived
            derived_rates::rate_Cr54_to_n_Cr53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr54_to_He4_Ti50_derived
            derived_rates::rate_Cr54_to_He4_Ti50_derived(rate_eval, tfactors, pf_cache);
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
            // Mn51_to_n_Mn50_derived
            derived_rates::rate_Mn51_to_n_Mn50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn51_to_p_Cr50_derived
            derived_rates::rate_Mn51_to_p_Cr50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn51_to_He4_V47_derived
            derived_rates::rate_Mn51_to_He4_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn52_to_n_Mn51_derived
            derived_rates::rate_Mn52_to_n_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn52_to_p_Cr51_derived
            derived_rates::rate_Mn52_to_p_Cr51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn52_to_He4_V48_derived
            derived_rates::rate_Mn52_to_He4_V48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn53_to_n_Mn52_derived
            derived_rates::rate_Mn53_to_n_Mn52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn53_to_p_Cr52_derived
            derived_rates::rate_Mn53_to_p_Cr52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn53_to_He4_V49_derived
            derived_rates::rate_Mn53_to_He4_V49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn54_to_n_Mn53_derived
            derived_rates::rate_Mn54_to_n_Mn53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn54_to_p_Cr53_derived
            derived_rates::rate_Mn54_to_p_Cr53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn54_to_He4_V50_derived
            derived_rates::rate_Mn54_to_He4_V50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn55_to_n_Mn54_derived
            derived_rates::rate_Mn55_to_n_Mn54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn55_to_p_Cr54_derived
            derived_rates::rate_Mn55_to_p_Cr54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn55_to_He4_V51_derived
            derived_rates::rate_Mn55_to_He4_V51_derived(rate_eval, tfactors, pf_cache);
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
            // Fe53_to_n_Fe52_derived
            derived_rates::rate_Fe53_to_n_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe53_to_p_Mn52_derived
            derived_rates::rate_Fe53_to_p_Mn52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe53_to_He4_Cr49_derived
            derived_rates::rate_Fe53_to_He4_Cr49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe54_to_n_Fe53_derived
            derived_rates::rate_Fe54_to_n_Fe53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe54_to_p_Mn53_derived
            derived_rates::rate_Fe54_to_p_Mn53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe54_to_He4_Cr50_derived
            derived_rates::rate_Fe54_to_He4_Cr50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe55_to_n_Fe54_derived
            derived_rates::rate_Fe55_to_n_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe55_to_p_Mn54_derived
            derived_rates::rate_Fe55_to_p_Mn54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe55_to_He4_Cr51_derived
            derived_rates::rate_Fe55_to_He4_Cr51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe56_to_n_Fe55_derived
            derived_rates::rate_Fe56_to_n_Fe55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe56_to_p_Mn55_derived
            derived_rates::rate_Fe56_to_p_Mn55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe56_to_He4_Cr52_derived
            derived_rates::rate_Fe56_to_He4_Cr52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe57_to_n_Fe56_derived
            derived_rates::rate_Fe57_to_n_Fe56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe57_to_He4_Cr53_derived
            derived_rates::rate_Fe57_to_He4_Cr53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe58_to_n_Fe57_derived
            derived_rates::rate_Fe58_to_n_Fe57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe58_to_He4_Cr54_derived
            derived_rates::rate_Fe58_to_He4_Cr54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co53_to_p_Fe52_derived
            derived_rates::rate_Co53_to_p_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co54_to_n_Co53_derived
            derived_rates::rate_Co54_to_n_Co53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co54_to_p_Fe53_derived
            derived_rates::rate_Co54_to_p_Fe53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co54_to_He4_Mn50_derived
            derived_rates::rate_Co54_to_He4_Mn50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_n_Co54_derived
            derived_rates::rate_Co55_to_n_Co54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_p_Fe54_derived
            derived_rates::rate_Co55_to_p_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_He4_Mn51_derived
            derived_rates::rate_Co55_to_He4_Mn51_derived(rate_eval, tfactors, pf_cache);
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
            // Co57_to_He4_Mn53_derived
            derived_rates::rate_Co57_to_He4_Mn53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co58_to_n_Co57_derived
            derived_rates::rate_Co58_to_n_Co57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co58_to_p_Fe57_derived
            derived_rates::rate_Co58_to_p_Fe57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co58_to_He4_Mn54_derived
            derived_rates::rate_Co58_to_He4_Mn54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co59_to_n_Co58_derived
            derived_rates::rate_Co59_to_n_Co58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co59_to_p_Fe58_derived
            derived_rates::rate_Co59_to_p_Fe58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co59_to_He4_Mn55_derived
            derived_rates::rate_Co59_to_He4_Mn55_derived(rate_eval, tfactors, pf_cache);
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
            // Ni59_to_n_Ni58_derived
            derived_rates::rate_Ni59_to_n_Ni58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni59_to_p_Co58_derived
            derived_rates::rate_Ni59_to_p_Co58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni59_to_He4_Fe55_derived
            derived_rates::rate_Ni59_to_He4_Fe55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni60_to_n_Ni59_derived
            derived_rates::rate_Ni60_to_n_Ni59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni60_to_p_Co59_derived
            derived_rates::rate_Ni60_to_p_Co59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni60_to_He4_Fe56_derived
            derived_rates::rate_Ni60_to_He4_Fe56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni61_to_n_Ni60_derived
            derived_rates::rate_Ni61_to_n_Ni60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni61_to_He4_Fe57_derived
            derived_rates::rate_Ni61_to_He4_Fe57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni62_to_n_Ni61_derived
            derived_rates::rate_Ni62_to_n_Ni61_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni62_to_He4_Fe58_derived
            derived_rates::rate_Ni62_to_He4_Fe58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni63_to_n_Ni62_derived
            derived_rates::rate_Ni63_to_n_Ni62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni64_to_n_Ni63_derived
            derived_rates::rate_Ni64_to_n_Ni63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu57_to_p_Ni56_derived
            derived_rates::rate_Cu57_to_p_Ni56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu57_to_He4_Co53_derived
            derived_rates::rate_Cu57_to_He4_Co53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu58_to_n_Cu57_derived
            derived_rates::rate_Cu58_to_n_Cu57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu58_to_p_Ni57_derived
            derived_rates::rate_Cu58_to_p_Ni57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu58_to_He4_Co54_derived
            derived_rates::rate_Cu58_to_He4_Co54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu59_to_n_Cu58_derived
            derived_rates::rate_Cu59_to_n_Cu58_derived(rate_eval, tfactors, pf_cache);
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
            // Cu60_to_n_Cu59_derived
            derived_rates::rate_Cu60_to_n_Cu59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu60_to_p_Ni59_derived
            derived_rates::rate_Cu60_to_p_Ni59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu60_to_He4_Co56_derived
            derived_rates::rate_Cu60_to_He4_Co56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu61_to_n_Cu60_derived
            derived_rates::rate_Cu61_to_n_Cu60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu61_to_p_Ni60_derived
            derived_rates::rate_Cu61_to_p_Ni60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu61_to_He4_Co57_derived
            derived_rates::rate_Cu61_to_He4_Co57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu62_to_n_Cu61_derived
            derived_rates::rate_Cu62_to_n_Cu61_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu62_to_p_Ni61_derived
            derived_rates::rate_Cu62_to_p_Ni61_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu62_to_He4_Co58_derived
            derived_rates::rate_Cu62_to_He4_Co58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu63_to_n_Cu62_derived
            derived_rates::rate_Cu63_to_n_Cu62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu63_to_p_Ni62_derived
            derived_rates::rate_Cu63_to_p_Ni62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu63_to_He4_Co59_derived
            derived_rates::rate_Cu63_to_He4_Co59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu64_to_n_Cu63_derived
            derived_rates::rate_Cu64_to_n_Cu63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu64_to_p_Ni63_derived
            derived_rates::rate_Cu64_to_p_Ni63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu65_to_n_Cu64_derived
            derived_rates::rate_Cu65_to_n_Cu64_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu65_to_p_Ni64_derived
            derived_rates::rate_Cu65_to_p_Ni64_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn59_to_p_Cu58_derived
            derived_rates::rate_Zn59_to_p_Cu58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn60_to_n_Zn59_derived
            derived_rates::rate_Zn60_to_n_Zn59_derived(rate_eval, tfactors, pf_cache);
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
            // Zn61_to_n_Zn60_derived
            derived_rates::rate_Zn61_to_n_Zn60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn61_to_p_Cu60_derived
            derived_rates::rate_Zn61_to_p_Cu60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn61_to_He4_Ni57_derived
            derived_rates::rate_Zn61_to_He4_Ni57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn62_to_n_Zn61_derived
            derived_rates::rate_Zn62_to_n_Zn61_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn62_to_p_Cu61_derived
            derived_rates::rate_Zn62_to_p_Cu61_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn62_to_He4_Ni58_derived
            derived_rates::rate_Zn62_to_He4_Ni58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn63_to_n_Zn62_derived
            derived_rates::rate_Zn63_to_n_Zn62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn63_to_p_Cu62_derived
            derived_rates::rate_Zn63_to_p_Cu62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn63_to_He4_Ni59_derived
            derived_rates::rate_Zn63_to_He4_Ni59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn64_to_n_Zn63_derived
            derived_rates::rate_Zn64_to_n_Zn63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn64_to_p_Cu63_derived
            derived_rates::rate_Zn64_to_p_Cu63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn64_to_He4_Ni60_derived
            derived_rates::rate_Zn64_to_He4_Ni60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn65_to_n_Zn64_derived
            derived_rates::rate_Zn65_to_n_Zn64_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn65_to_p_Cu64_derived
            derived_rates::rate_Zn65_to_p_Cu64_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn65_to_He4_Ni61_derived
            derived_rates::rate_Zn65_to_He4_Ni61_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn66_to_n_Zn65_derived
            derived_rates::rate_Zn66_to_n_Zn65_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn66_to_p_Cu65_derived
            derived_rates::rate_Zn66_to_p_Cu65_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn66_to_He4_Ni62_derived
            derived_rates::rate_Zn66_to_He4_Ni62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga62_to_p_Zn61_derived
            derived_rates::rate_Ga62_to_p_Zn61_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga62_to_He4_Cu58_derived
            derived_rates::rate_Ga62_to_He4_Cu58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga63_to_n_Ga62_derived
            derived_rates::rate_Ga63_to_n_Ga62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga63_to_p_Zn62_derived
            derived_rates::rate_Ga63_to_p_Zn62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga63_to_He4_Cu59_derived
            derived_rates::rate_Ga63_to_He4_Cu59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga64_to_n_Ga63_derived
            derived_rates::rate_Ga64_to_n_Ga63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga64_to_p_Zn63_derived
            derived_rates::rate_Ga64_to_p_Zn63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ga64_to_He4_Cu60_derived
            derived_rates::rate_Ga64_to_He4_Cu60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ge63_to_p_Ga62_derived
            derived_rates::rate_Ge63_to_p_Ga62_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ge63_to_He4_Zn59_derived
            derived_rates::rate_Ge63_to_He4_Zn59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ge64_to_n_Ge63_derived
            derived_rates::rate_Ge64_to_n_Ge63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ge64_to_p_Ga63_derived
            derived_rates::rate_Ge64_to_p_Ga63_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ge64_to_He4_Zn60_derived
            derived_rates::rate_Ge64_to_He4_Zn60_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Li6_to_n_p_He4_derived
            derived_rates::rate_Li6_to_n_p_He4_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Be9_to_n_He4_He4_derived
            derived_rates::rate_Be9_to_n_He4_He4_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C12_to_3He4_derived
            derived_rates::rate_C12_to_3He4_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // He3_n_to_d_d_derived
            derived_rates::rate_He3_n_to_d_d_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // He3_He4_to_p_Li6_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He3_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He3_He4);
            }
#endif
            derived_rates::rate_He3_He4_to_p_Li6_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // He4_He4_to_n_Be7_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4);
            }
#endif
            derived_rates::rate_He4_He4_to_n_Be7_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_He4_to_p_Li7_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4);
            }
#endif
            derived_rates::rate_He4_He4_to_p_Li7_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Li6_He4_to_p_Be9_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Li6);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Li6);
            }
#endif
            derived_rates::rate_Li6_He4_to_p_Be9_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Li7_p_to_n_Be7_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Li7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Li7);
            }
#endif
            derived_rates::rate_Li7_p_to_n_Be7_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Li7_p_to_d_Li6_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Li7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Li7);
            }
#endif
            derived_rates::rate_Li7_p_to_d_Li6_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Li7_He4_to_n_B10_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Li7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Li7);
            }
#endif
            derived_rates::rate_Li7_He4_to_n_B10_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Be7_n_to_d_Li6_derived
            derived_rates::rate_Be7_n_to_d_Li6_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Be7_He4_to_p_B10_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Be7);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Be7);
            }
#endif
            derived_rates::rate_Be7_He4_to_p_B10_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // C12_n_to_He4_Be9_derived
            derived_rates::rate_C12_n_to_He4_Be9_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C12_He4_to_n_O15_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C12);
            }
#endif
            derived_rates::rate_C12_He4_to_n_O15_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // C12_C12_to_n_Mg23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            derived_rates::rate_C12_C12_to_n_Mg23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // C13_p_to_He4_B10_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_C13);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C13);
            }
#endif
            derived_rates::rate_C13_p_to_He4_B10_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // C14_p_to_n_N14_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_C14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C14);
            }
#endif
            derived_rates::rate_C14_p_to_n_N14_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // C14_p_to_He4_B11_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_C14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C14);
            }
#endif
            derived_rates::rate_C14_p_to_He4_B11_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // C14_He4_to_n_O17_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C14);
            }
#endif
            derived_rates::rate_C14_He4_to_n_O17_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // N13_n_to_p_C13_derived
            derived_rates::rate_N13_n_to_p_C13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N13_n_to_He4_B10_derived
            derived_rates::rate_N13_n_to_He4_B10_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N14_n_to_d_C13_derived
            derived_rates::rate_N14_n_to_d_C13_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N14_n_to_He4_B11_derived
            derived_rates::rate_N14_n_to_He4_B11_derived(rate_eval, tfactors, pf_cache);
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
            // N15_n_to_d_C14_derived
            derived_rates::rate_N15_n_to_d_C14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // N15_He4_to_n_F18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N15);
            }
#endif
            derived_rates::rate_N15_He4_to_n_F18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // O14_n_to_p_N14_derived
            derived_rates::rate_O14_n_to_p_N14_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O15_n_to_p_N15_derived
            derived_rates::rate_O15_n_to_p_N15_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // O15_He4_to_n_Ne18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O15);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O15);
            }
#endif
            derived_rates::rate_O15_He4_to_n_Ne18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // O16_n_to_He4_C13_derived
            derived_rates::rate_O16_n_to_He4_C13_derived(rate_eval, tfactors, pf_cache);
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
            // O16_He4_to_n_Ne19_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O16);
            }
#endif
            derived_rates::rate_O16_He4_to_n_Ne19_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // O17_p_to_n_F17_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O17);
            }
#endif
            derived_rates::rate_O17_p_to_n_F17_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // O18_p_to_n_F18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O18);
            }
#endif
            derived_rates::rate_O18_p_to_n_F18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // F17_n_to_He4_N14_derived
            derived_rates::rate_F17_n_to_He4_N14_derived(rate_eval, tfactors, pf_cache);
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
            // F17_He4_to_p_Ne20_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F17);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F17);
            }
#endif
            derived_rates::rate_F17_He4_to_p_Ne20_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // F18_p_to_n_Ne18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_F18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F18);
            }
#endif
            derived_rates::rate_F18_p_to_n_Ne18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // F18_He4_to_n_Na21_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F18);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F18);
            }
#endif
            derived_rates::rate_F18_He4_to_n_Na21_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // F19_He4_to_n_Na22_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_F19);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F19);
            }
#endif
            derived_rates::rate_F19_He4_to_n_Na22_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ne19_n_to_p_F19_derived
            derived_rates::rate_Ne19_n_to_p_F19_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_n_to_He4_O17_derived
            derived_rates::rate_Ne20_n_to_He4_O17_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_He4_to_n_Mg23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            derived_rates::rate_Ne20_He4_to_n_Mg23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Ne21_n_to_He4_O18_derived
            derived_rates::rate_Ne21_n_to_He4_O18_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne21_p_to_n_Na21_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne21);
            }
#endif
            derived_rates::rate_Ne21_p_to_n_Na21_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Ne22_p_to_n_Na22_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ne22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne22);
            }
#endif
            derived_rates::rate_Ne22_p_to_n_Na22_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Na21_p_to_He4_Ne18_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na21);
            }
#endif
            derived_rates::rate_Na21_p_to_He4_Ne18_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Na21_He4_to_p_Mg24_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na21);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na21);
            }
#endif
            derived_rates::rate_Na21_He4_to_p_Mg24_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Na22_p_to_He4_Ne19_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na22);
            }
#endif
            derived_rates::rate_Na22_p_to_He4_Ne19_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Na22_He4_to_n_Al25_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na22);
            }
#endif
            derived_rates::rate_Na22_He4_to_n_Al25_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Na23_He4_to_n_Al26_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na23);
            }
#endif
            derived_rates::rate_Na23_He4_to_n_Al26_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg23_n_to_p_Na23_derived
            derived_rates::rate_Mg23_n_to_p_Na23_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg24_n_to_He4_Ne21_derived
            derived_rates::rate_Mg24_n_to_He4_Ne21_derived(rate_eval, tfactors, pf_cache);
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
            // Mg25_n_to_He4_Ne22_derived
            derived_rates::rate_Mg25_n_to_He4_Ne22_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg25_p_to_n_Al25_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg25);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg25);
            }
#endif
            derived_rates::rate_Mg25_p_to_n_Al25_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Mg26_p_to_n_Al26_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg26);
            }
#endif
            derived_rates::rate_Mg26_p_to_n_Al26_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg26_p_to_He4_Na23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg26);
            }
#endif
            derived_rates::rate_Mg26_p_to_He4_Na23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Al26_p_to_He4_Mg23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al26);
            }
#endif
            derived_rates::rate_Al26_p_to_He4_Mg23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Al26_He4_to_n_P29_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al26);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al26);
            }
#endif
            derived_rates::rate_Al26_He4_to_n_P29_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Si28_n_to_He4_Mg25_derived
            derived_rates::rate_Si28_n_to_He4_Mg25_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_p_to_He4_Al25_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si28);
            }
#endif
            derived_rates::rate_Si28_p_to_He4_Al25_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Si29_n_to_He4_Mg26_derived
            derived_rates::rate_Si29_n_to_He4_Mg26_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si29_p_to_n_P29_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si29);
            }
#endif
            derived_rates::rate_Si29_p_to_n_P29_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Si29_He4_to_n_S32_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si29);
            }
#endif
            derived_rates::rate_Si29_He4_to_n_S32_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si29_He4_to_p_P32_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si29);
            }
#endif
            derived_rates::rate_Si29_He4_to_p_P32_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si30_p_to_n_P30_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si30);
            }
#endif
            derived_rates::rate_Si30_p_to_n_P30_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si30_p_to_He4_Al27_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si30);
            }
#endif
            derived_rates::rate_Si30_p_to_He4_Al27_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si30_He4_to_n_S33_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si30);
            }
#endif
            derived_rates::rate_Si30_He4_to_n_S33_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si30_He4_to_p_P33_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si30);
            }
#endif
            derived_rates::rate_Si30_He4_to_p_P33_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si32_p_to_n_P32_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si32);
            }
#endif
            derived_rates::rate_Si32_p_to_n_P32_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si32_He4_to_n_S35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si32);
            }
#endif
            derived_rates::rate_Si32_He4_to_n_S35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P30_n_to_He4_Al27_derived
            derived_rates::rate_P30_n_to_He4_Al27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P30_He4_to_n_Cl33_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P30);
            }
#endif
            derived_rates::rate_P30_He4_to_n_Cl33_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P31_n_to_p_Si31_derived
            derived_rates::rate_P31_n_to_p_Si31_derived(rate_eval, tfactors, pf_cache);
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
            // P31_He4_to_n_Cl34_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P31);
            }
#endif
            derived_rates::rate_P31_He4_to_n_Cl34_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P32_He4_to_n_Cl35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P32);
            }
#endif
            derived_rates::rate_P32_He4_to_n_Cl35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P32_He4_to_p_S35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P32);
            }
#endif
            derived_rates::rate_P32_He4_to_p_S35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P33_p_to_n_S33_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P33);
            }
#endif
            derived_rates::rate_P33_p_to_n_S33_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P33_He4_to_n_Cl36_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P33);
            }
#endif
            derived_rates::rate_P33_He4_to_n_Cl36_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P33_He4_to_p_S36_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P33);
            }
#endif
            derived_rates::rate_P33_He4_to_p_S36_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S32_n_to_p_P32_derived
            derived_rates::rate_S32_n_to_p_P32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_p_to_He4_P29_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S32);
            }
#endif
            derived_rates::rate_S32_p_to_He4_P29_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // S33_p_to_n_Cl33_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S33);
            }
#endif
            derived_rates::rate_S33_p_to_n_Cl33_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // S33_He4_to_n_Ar36_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S33);
            }
#endif
            derived_rates::rate_S33_He4_to_n_Ar36_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S33_He4_to_p_Cl36_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S33);
            }
#endif
            derived_rates::rate_S33_He4_to_p_Cl36_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S34_n_to_He4_Si31_derived
            derived_rates::rate_S34_n_to_He4_Si31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S34_p_to_n_Cl34_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S34);
            }
#endif
            derived_rates::rate_S34_p_to_n_Cl34_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S34_p_to_He4_P31_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S34);
            }
#endif
            derived_rates::rate_S34_p_to_He4_P31_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S35_p_to_n_Cl35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S35);
            }
#endif
            derived_rates::rate_S35_p_to_n_Cl35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S36_p_to_n_Cl36_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S36);
            }
#endif
            derived_rates::rate_S36_p_to_n_Cl36_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S36_He4_to_n_Ar39_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S36);
            }
#endif
            derived_rates::rate_S36_He4_to_n_Ar39_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl34_He4_to_n_K37_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl34);
            }
#endif
            derived_rates::rate_Cl34_He4_to_n_K37_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl35_He4_to_n_K38_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl35);
            }
#endif
            derived_rates::rate_Cl35_He4_to_n_K38_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl36_p_to_n_Ar36_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl36);
            }
#endif
            derived_rates::rate_Cl36_p_to_n_Ar36_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl36_He4_to_n_K39_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl36);
            }
#endif
            derived_rates::rate_Cl36_He4_to_n_K39_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl36_He4_to_p_Ar39_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl36);
            }
#endif
            derived_rates::rate_Cl36_He4_to_p_Ar39_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl37_p_to_n_Ar37_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl37);
            }
#endif
            derived_rates::rate_Cl37_p_to_n_Ar37_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl37_p_to_He4_S34_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl37);
            }
#endif
            derived_rates::rate_Cl37_p_to_He4_S34_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl37_He4_to_n_K40_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl37);
            }
#endif
            derived_rates::rate_Cl37_He4_to_n_K40_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl37_He4_to_p_Ar40_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl37);
            }
#endif
            derived_rates::rate_Cl37_He4_to_p_Ar40_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar36_p_to_He4_Cl33_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar36);
            }
#endif
            derived_rates::rate_Ar36_p_to_He4_Cl33_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Ar37_n_to_He4_S34_derived
            derived_rates::rate_Ar37_n_to_He4_S34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar37_p_to_n_K37_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar37);
            }
#endif
            derived_rates::rate_Ar37_p_to_n_K37_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Ar37_He4_to_n_Ca40_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar37);
            }
#endif
            derived_rates::rate_Ar37_He4_to_n_Ca40_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar37_He4_to_p_K40_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar37);
            }
#endif
            derived_rates::rate_Ar37_He4_to_p_K40_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar38_n_to_He4_S35_derived
            derived_rates::rate_Ar38_n_to_He4_S35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar38_p_to_n_K38_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar38);
            }
#endif
            derived_rates::rate_Ar38_p_to_n_K38_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar38_p_to_He4_Cl35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar38);
            }
#endif
            derived_rates::rate_Ar38_p_to_He4_Cl35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar38_He4_to_n_Ca41_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar38);
            }
#endif
            derived_rates::rate_Ar38_He4_to_n_Ca41_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar38_He4_to_p_K41_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar38);
            }
#endif
            derived_rates::rate_Ar38_He4_to_p_K41_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar39_p_to_n_K39_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar39);
            }
#endif
            derived_rates::rate_Ar39_p_to_n_K39_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar39_He4_to_n_Ca42_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar39);
            }
#endif
            derived_rates::rate_Ar39_He4_to_n_Ca42_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar40_p_to_n_K40_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar40);
            }
#endif
            derived_rates::rate_Ar40_p_to_n_K40_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar40_He4_to_n_Ca43_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar40);
            }
#endif
            derived_rates::rate_Ar40_He4_to_n_Ca43_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // K39_He4_to_p_Ca42_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K39);
            }
#endif
            derived_rates::rate_K39_He4_to_p_Ca42_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // K40_He4_to_n_Sc43_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K40);
            }
#endif
            derived_rates::rate_K40_He4_to_n_Sc43_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca40_n_to_p_K40_derived
            derived_rates::rate_Ca40_n_to_p_K40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_p_to_He4_K37_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca40);
            }
#endif
            derived_rates::rate_Ca40_p_to_He4_K37_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Ca41_n_to_p_K41_derived
            derived_rates::rate_Ca41_n_to_p_K41_derived(rate_eval, tfactors, pf_cache);
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
            // Ca41_He4_to_n_Ti44_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca41);
            }
#endif
            derived_rates::rate_Ca41_He4_to_n_Ti44_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca41_He4_to_p_Sc44_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca41);
            }
#endif
            derived_rates::rate_Ca41_He4_to_p_Sc44_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca42_He4_to_n_Ti45_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca42);
            }
#endif
            derived_rates::rate_Ca42_He4_to_n_Ti45_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca42_He4_to_p_Sc45_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca42);
            }
#endif
            derived_rates::rate_Ca42_He4_to_p_Sc45_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca43_p_to_n_Sc43_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca43);
            }
#endif
            derived_rates::rate_Ca43_p_to_n_Sc43_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca43_p_to_He4_K40_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca43);
            }
#endif
            derived_rates::rate_Ca43_p_to_He4_K40_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca43_He4_to_p_Sc46_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca43);
            }
#endif
            derived_rates::rate_Ca43_He4_to_p_Sc46_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca44_p_to_n_Sc44_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca44);
            }
#endif
            derived_rates::rate_Ca44_p_to_n_Sc44_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca44_p_to_He4_K41_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca44);
            }
#endif
            derived_rates::rate_Ca44_p_to_He4_K41_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca44_He4_to_n_Ti47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca44);
            }
#endif
            derived_rates::rate_Ca44_He4_to_n_Ti47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca44_He4_to_p_Sc47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca44);
            }
#endif
            derived_rates::rate_Ca44_He4_to_p_Sc47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca45_p_to_n_Sc45_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca45);
            }
#endif
            derived_rates::rate_Ca45_p_to_n_Sc45_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca45_He4_to_p_Sc48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca45);
            }
#endif
            derived_rates::rate_Ca45_He4_to_p_Sc48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca46_p_to_n_Sc46_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca46);
            }
#endif
            derived_rates::rate_Ca46_p_to_n_Sc46_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca46_He4_to_n_Ti49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca46);
            }
#endif
            derived_rates::rate_Ca46_He4_to_n_Ti49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca46_He4_to_p_Sc49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca46);
            }
#endif
            derived_rates::rate_Ca46_He4_to_p_Sc49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca48_p_to_n_Sc48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca48);
            }
#endif
            derived_rates::rate_Ca48_p_to_n_Sc48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca48_He4_to_n_Ti51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca48);
            }
#endif
            derived_rates::rate_Ca48_He4_to_n_Ti51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc43_He4_to_n_V46_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc43);
            }
#endif
            derived_rates::rate_Sc43_He4_to_n_V46_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc44_n_to_He4_K41_derived
            derived_rates::rate_Sc44_n_to_He4_K41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc44_p_to_n_Ti44_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc44);
            }
#endif
            derived_rates::rate_Sc44_p_to_n_Ti44_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc44_He4_to_n_V47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc44);
            }
#endif
            derived_rates::rate_Sc44_He4_to_n_V47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc45_p_to_n_Ti45_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc45);
            }
#endif
            derived_rates::rate_Sc45_p_to_n_Ti45_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc45_He4_to_n_V48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc45);
            }
#endif
            derived_rates::rate_Sc45_He4_to_n_V48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc47_n_to_p_Ca47_derived
            derived_rates::rate_Sc47_n_to_p_Ca47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc47_p_to_n_Ti47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc47);
            }
#endif
            derived_rates::rate_Sc47_p_to_n_Ti47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc47_He4_to_n_V50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc47);
            }
#endif
            derived_rates::rate_Sc47_He4_to_n_V50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc49_He4_to_n_V52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc49);
            }
#endif
            derived_rates::rate_Sc49_He4_to_n_V52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti45_He4_to_n_Cr48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti45);
            }
#endif
            derived_rates::rate_Ti45_He4_to_n_Cr48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti46_n_to_p_Sc46_derived
            derived_rates::rate_Ti46_n_to_p_Sc46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti46_n_to_He4_Ca43_derived
            derived_rates::rate_Ti46_n_to_He4_Ca43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti46_p_to_n_V46_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti46);
            }
#endif
            derived_rates::rate_Ti46_p_to_n_V46_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti46_p_to_He4_Sc43_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti46);
            }
#endif
            derived_rates::rate_Ti46_p_to_He4_Sc43_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti46_He4_to_n_Cr49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti46);
            }
#endif
            derived_rates::rate_Ti46_He4_to_n_Cr49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti46_He4_to_p_V49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti46);
            }
#endif
            derived_rates::rate_Ti46_He4_to_p_V49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti47_p_to_n_V47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti47);
            }
#endif
            derived_rates::rate_Ti47_p_to_n_V47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti47_p_to_He4_Sc44_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti47);
            }
#endif
            derived_rates::rate_Ti47_p_to_He4_Sc44_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti47_He4_to_n_Cr50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti47);
            }
#endif
            derived_rates::rate_Ti47_He4_to_n_Cr50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti47_He4_to_p_V50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti47);
            }
#endif
            derived_rates::rate_Ti47_He4_to_p_V50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti48_n_to_p_Sc48_derived
            derived_rates::rate_Ti48_n_to_p_Sc48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti48_n_to_He4_Ca45_derived
            derived_rates::rate_Ti48_n_to_He4_Ca45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti48_p_to_n_V48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti48);
            }
#endif
            derived_rates::rate_Ti48_p_to_n_V48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti48_p_to_He4_Sc45_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti48);
            }
#endif
            derived_rates::rate_Ti48_p_to_He4_Sc45_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti48_He4_to_n_Cr51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti48);
            }
#endif
            derived_rates::rate_Ti48_He4_to_n_Cr51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti48_He4_to_p_V51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti48);
            }
#endif
            derived_rates::rate_Ti48_He4_to_p_V51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti49_n_to_p_Sc49_derived
            derived_rates::rate_Ti49_n_to_p_Sc49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti49_p_to_n_V49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti49);
            }
#endif
            derived_rates::rate_Ti49_p_to_n_V49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti49_p_to_He4_Sc46_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti49);
            }
#endif
            derived_rates::rate_Ti49_p_to_He4_Sc46_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti49_He4_to_p_V52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti49);
            }
#endif
            derived_rates::rate_Ti49_He4_to_p_V52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti50_n_to_He4_Ca47_derived
            derived_rates::rate_Ti50_n_to_He4_Ca47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti50_p_to_n_V50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti50);
            }
#endif
            derived_rates::rate_Ti50_p_to_n_V50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti50_p_to_He4_Sc47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti50);
            }
#endif
            derived_rates::rate_Ti50_p_to_He4_Sc47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti50_He4_to_n_Cr53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti50);
            }
#endif
            derived_rates::rate_Ti50_He4_to_n_Cr53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti51_p_to_He4_Sc48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti51);
            }
#endif
            derived_rates::rate_Ti51_p_to_He4_Sc48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // V47_He4_to_n_Mn50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V47);
            }
#endif
            derived_rates::rate_V47_He4_to_n_Mn50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V48_p_to_n_Cr48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V48);
            }
#endif
            derived_rates::rate_V48_p_to_n_Cr48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V48_p_to_He4_Ti45_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V48);
            }
#endif
            derived_rates::rate_V48_p_to_He4_Ti45_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V48_He4_to_n_Mn51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V48);
            }
#endif
            derived_rates::rate_V48_He4_to_n_Mn51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V49_n_to_He4_Sc46_derived
            derived_rates::rate_V49_n_to_He4_Sc46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V49_p_to_n_Cr49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V49);
            }
#endif
            derived_rates::rate_V49_p_to_n_Cr49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V49_He4_to_n_Mn52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V49);
            }
#endif
            derived_rates::rate_V49_He4_to_n_Mn52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V50_He4_to_n_Mn53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V50);
            }
#endif
            derived_rates::rate_V50_He4_to_n_Mn53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V51_n_to_p_Ti51_derived
            derived_rates::rate_V51_n_to_p_Ti51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V51_n_to_He4_Sc48_derived
            derived_rates::rate_V51_n_to_He4_Sc48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V51_p_to_n_Cr51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V51);
            }
#endif
            derived_rates::rate_V51_p_to_n_Cr51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V51_He4_to_n_Mn54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V51);
            }
#endif
            derived_rates::rate_V51_He4_to_n_Mn54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V51_He4_to_p_Cr54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V51);
            }
#endif
            derived_rates::rate_V51_He4_to_p_Cr54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Cr49_He4_to_n_Fe52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr49);
            }
#endif
            derived_rates::rate_Cr49_He4_to_n_Fe52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr50_n_to_p_V50_derived
            derived_rates::rate_Cr50_n_to_p_V50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr50_p_to_n_Mn50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr50);
            }
#endif
            derived_rates::rate_Cr50_p_to_n_Mn50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr50_p_to_He4_V47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr50);
            }
#endif
            derived_rates::rate_Cr50_p_to_He4_V47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr50_He4_to_n_Fe53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr50);
            }
#endif
            derived_rates::rate_Cr50_He4_to_n_Fe53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr50_He4_to_p_Mn53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr50);
            }
#endif
            derived_rates::rate_Cr50_He4_to_p_Mn53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr51_p_to_n_Mn51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr51);
            }
#endif
            derived_rates::rate_Cr51_p_to_n_Mn51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr51_p_to_He4_V48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr51);
            }
#endif
            derived_rates::rate_Cr51_p_to_He4_V48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr51_He4_to_n_Fe54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr51);
            }
#endif
            derived_rates::rate_Cr51_He4_to_n_Fe54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr51_He4_to_p_Mn54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr51);
            }
#endif
            derived_rates::rate_Cr51_He4_to_p_Mn54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr52_n_to_p_V52_derived
            derived_rates::rate_Cr52_n_to_p_V52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr52_n_to_He4_Ti49_derived
            derived_rates::rate_Cr52_n_to_He4_Ti49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr52_p_to_n_Mn52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr52);
            }
#endif
            derived_rates::rate_Cr52_p_to_n_Mn52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr52_p_to_He4_V49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr52);
            }
#endif
            derived_rates::rate_Cr52_p_to_He4_V49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr52_He4_to_n_Fe55_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr52);
            }
#endif
            derived_rates::rate_Cr52_He4_to_n_Fe55_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr52_He4_to_p_Mn55_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr52);
            }
#endif
            derived_rates::rate_Cr52_He4_to_p_Mn55_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr53_p_to_n_Mn53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr53);
            }
#endif
            derived_rates::rate_Cr53_p_to_n_Mn53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr53_p_to_He4_V50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr53);
            }
#endif
            derived_rates::rate_Cr53_p_to_He4_V50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr53_He4_to_n_Fe56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr53);
            }
#endif
            derived_rates::rate_Cr53_He4_to_n_Fe56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr54_n_to_He4_Ti51_derived
            derived_rates::rate_Cr54_n_to_He4_Ti51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr54_p_to_n_Mn54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr54);
            }
#endif
            derived_rates::rate_Cr54_p_to_n_Mn54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr54_He4_to_n_Fe57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr54);
            }
#endif
            derived_rates::rate_Cr54_He4_to_n_Fe57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn50_He4_to_n_Co53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn50);
            }
#endif
            derived_rates::rate_Mn50_He4_to_n_Co53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Mn51_He4_to_n_Co54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn51);
            }
#endif
            derived_rates::rate_Mn51_He4_to_n_Co54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Mn53_p_to_n_Fe53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn53);
            }
#endif
            derived_rates::rate_Mn53_p_to_n_Fe53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn53_He4_to_n_Co56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn53);
            }
#endif
            derived_rates::rate_Mn53_He4_to_n_Co56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn54_p_to_n_Fe54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn54);
            }
#endif
            derived_rates::rate_Mn54_p_to_n_Fe54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn54_He4_to_n_Co57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn54);
            }
#endif
            derived_rates::rate_Mn54_He4_to_n_Co57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn54_He4_to_p_Fe57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn54);
            }
#endif
            derived_rates::rate_Mn54_He4_to_p_Fe57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn55_n_to_He4_V52_derived
            derived_rates::rate_Mn55_n_to_He4_V52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn55_p_to_n_Fe55_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn55);
            }
#endif
            derived_rates::rate_Mn55_p_to_n_Fe55_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn55_He4_to_n_Co58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn55);
            }
#endif
            derived_rates::rate_Mn55_He4_to_n_Co58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn55_He4_to_p_Fe58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn55);
            }
#endif
            derived_rates::rate_Mn55_He4_to_p_Fe58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe53_p_to_n_Co53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe53);
            }
#endif
            derived_rates::rate_Fe53_p_to_n_Co53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Fe54_p_to_n_Co54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
            }
#endif
            derived_rates::rate_Fe54_p_to_n_Co54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Fe55_He4_to_p_Co58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe55);
            }
#endif
            derived_rates::rate_Fe55_He4_to_p_Co58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Fe56_p_to_He4_Mn53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            derived_rates::rate_Fe56_p_to_He4_Mn53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe56_He4_to_n_Ni59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe56);
            }
#endif
            derived_rates::rate_Fe56_He4_to_n_Ni59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe56_He4_to_p_Co59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe56);
            }
#endif
            derived_rates::rate_Fe56_He4_to_p_Co59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe57_p_to_n_Co57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe57);
            }
#endif
            derived_rates::rate_Fe57_p_to_n_Co57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe57_He4_to_n_Ni60_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe57);
            }
#endif
            derived_rates::rate_Fe57_He4_to_n_Ni60_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe58_p_to_n_Co58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe58);
            }
#endif
            derived_rates::rate_Fe58_p_to_n_Co58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe58_He4_to_n_Ni61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe58);
            }
#endif
            derived_rates::rate_Fe58_He4_to_n_Ni61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co54_He4_to_n_Cu57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co54);
            }
#endif
            derived_rates::rate_Co54_He4_to_n_Cu57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Co55_He4_to_n_Cu58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co55);
            }
#endif
            derived_rates::rate_Co55_He4_to_n_Cu58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Co57_He4_to_n_Cu60_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co57);
            }
#endif
            derived_rates::rate_Co57_He4_to_n_Cu60_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co58_p_to_n_Ni58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co58);
            }
#endif
            derived_rates::rate_Co58_p_to_n_Ni58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co58_He4_to_n_Cu61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co58);
            }
#endif
            derived_rates::rate_Co58_He4_to_n_Cu61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co58_He4_to_p_Ni61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co58);
            }
#endif
            derived_rates::rate_Co58_He4_to_p_Ni61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co59_p_to_n_Ni59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co59);
            }
#endif
            derived_rates::rate_Co59_p_to_n_Ni59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co59_He4_to_n_Cu62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co59);
            }
#endif
            derived_rates::rate_Co59_He4_to_n_Cu62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co59_He4_to_p_Ni62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co59);
            }
#endif
            derived_rates::rate_Co59_He4_to_p_Ni62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni56_p_to_He4_Co53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni56);
            }
#endif
            derived_rates::rate_Ni56_p_to_He4_Co53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni56_He4_to_n_Zn59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
            }
#endif
            derived_rates::rate_Ni56_He4_to_n_Zn59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Ni57_p_to_n_Cu57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni57);
            }
#endif
            derived_rates::rate_Ni57_p_to_n_Cu57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni57_p_to_He4_Co54_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni57);
            }
#endif
            derived_rates::rate_Ni57_p_to_He4_Co54_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Ni57_He4_to_p_Cu60_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni57);
            }
#endif
            derived_rates::rate_Ni57_He4_to_p_Cu60_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni58_p_to_n_Cu58_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
            }
#endif
            derived_rates::rate_Ni58_p_to_n_Cu58_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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

        {
            // Ni58_He4_to_n_Zn61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni58);
            }
#endif
            derived_rates::rate_Ni58_He4_to_n_Zn61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni58_He4_to_p_Cu61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni58);
            }
#endif
            derived_rates::rate_Ni58_He4_to_p_Cu61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni59_p_to_n_Cu59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni59);
            }
#endif
            derived_rates::rate_Ni59_p_to_n_Cu59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni59_p_to_He4_Co56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni59);
            }
#endif
            derived_rates::rate_Ni59_p_to_He4_Co56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni59_He4_to_n_Zn62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni59);
            }
#endif
            derived_rates::rate_Ni59_He4_to_n_Zn62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni59_He4_to_p_Cu62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni59);
            }
#endif
            derived_rates::rate_Ni59_He4_to_p_Cu62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni60_p_to_n_Cu60_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni60);
            }
#endif
            derived_rates::rate_Ni60_p_to_n_Cu60_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni60_p_to_He4_Co57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni60);
            }
#endif
            derived_rates::rate_Ni60_p_to_He4_Co57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni60_He4_to_n_Zn63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni60);
            }
#endif
            derived_rates::rate_Ni60_He4_to_n_Zn63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni60_He4_to_p_Cu63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni60);
            }
#endif
            derived_rates::rate_Ni60_He4_to_p_Cu63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni61_p_to_n_Cu61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni61);
            }
#endif
            derived_rates::rate_Ni61_p_to_n_Cu61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni61_He4_to_n_Zn64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni61);
            }
#endif
            derived_rates::rate_Ni61_He4_to_n_Zn64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni61_He4_to_p_Cu64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni61);
            }
#endif
            derived_rates::rate_Ni61_He4_to_p_Cu64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni62_p_to_n_Cu62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni62);
            }
#endif
            derived_rates::rate_Ni62_p_to_n_Cu62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni62_He4_to_n_Zn65_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni62);
            }
#endif
            derived_rates::rate_Ni62_He4_to_n_Zn65_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni62_He4_to_p_Cu65_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni62);
            }
#endif
            derived_rates::rate_Ni62_He4_to_p_Cu65_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni63_p_to_n_Cu63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni63);
            }
#endif
            derived_rates::rate_Ni63_p_to_n_Cu63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni63_He4_to_n_Zn66_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni63);
            }
#endif
            derived_rates::rate_Ni63_He4_to_n_Zn66_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni64_p_to_n_Cu64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni64);
            }
#endif
            derived_rates::rate_Ni64_p_to_n_Cu64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu58_He4_to_p_Zn61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu58);
            }
#endif
            derived_rates::rate_Cu58_He4_to_p_Zn61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu59_p_to_n_Zn59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu59);
            }
#endif
            derived_rates::rate_Cu59_p_to_n_Zn59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu59_He4_to_n_Ga62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu59);
            }
#endif
            derived_rates::rate_Cu59_He4_to_n_Ga62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu59_He4_to_p_Zn62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu59);
            }
#endif
            derived_rates::rate_Cu59_He4_to_p_Zn62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu60_p_to_n_Zn60_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu60);
            }
#endif
            derived_rates::rate_Cu60_p_to_n_Zn60_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu60_He4_to_n_Ga63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu60);
            }
#endif
            derived_rates::rate_Cu60_He4_to_n_Ga63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu60_He4_to_p_Zn63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu60);
            }
#endif
            derived_rates::rate_Cu60_He4_to_p_Zn63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu61_p_to_n_Zn61_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu61);
            }
#endif
            derived_rates::rate_Cu61_p_to_n_Zn61_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu61_He4_to_n_Ga64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu61);
            }
#endif
            derived_rates::rate_Cu61_He4_to_n_Ga64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu61_He4_to_p_Zn64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu61);
            }
#endif
            derived_rates::rate_Cu61_He4_to_p_Zn64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu62_p_to_n_Zn62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu62);
            }
#endif
            derived_rates::rate_Cu62_p_to_n_Zn62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu62_He4_to_p_Zn65_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu62);
            }
#endif
            derived_rates::rate_Cu62_He4_to_p_Zn65_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu63_p_to_n_Zn63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu63);
            }
#endif
            derived_rates::rate_Cu63_p_to_n_Zn63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu63_He4_to_p_Zn66_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cu63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cu63);
            }
#endif
            derived_rates::rate_Cu63_He4_to_p_Zn66_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu64_p_to_n_Zn64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu64);
            }
#endif
            derived_rates::rate_Cu64_p_to_n_Zn64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cu65_p_to_n_Zn65_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu65);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu65);
            }
#endif
            derived_rates::rate_Cu65_p_to_n_Zn65_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn59_He4_to_p_Ga62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Zn59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Zn59);
            }
#endif
            derived_rates::rate_Zn59_He4_to_p_Ga62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn60_p_to_He4_Cu57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn60);
            }
#endif
            derived_rates::rate_Zn60_p_to_He4_Cu57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn60_He4_to_n_Ge63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Zn60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Zn60);
            }
#endif
            derived_rates::rate_Zn60_He4_to_n_Ge63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn60_He4_to_p_Ga63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Zn60);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Zn60);
            }
#endif
            derived_rates::rate_Zn60_He4_to_p_Ga63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn61_He4_to_n_Ge64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Zn61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Zn61);
            }
#endif
            derived_rates::rate_Zn61_He4_to_n_Ge64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn61_He4_to_p_Ga64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Zn61);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Zn61);
            }
#endif
            derived_rates::rate_Zn61_He4_to_p_Ga64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn62_p_to_n_Ga62_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn62);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn62);
            }
#endif
            derived_rates::rate_Zn62_p_to_n_Ga62_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn63_p_to_n_Ga63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn63);
            }
#endif
            derived_rates::rate_Zn63_p_to_n_Ga63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Zn64_p_to_n_Ga64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Zn64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Zn64);
            }
#endif
            derived_rates::rate_Zn64_p_to_n_Ga64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ga63_p_to_n_Ge63_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ga63);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ga63);
            }
#endif
            derived_rates::rate_Ga63_p_to_n_Ge63_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ga64_p_to_n_Ge64_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ga64);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ga64);
            }
#endif
            derived_rates::rate_Ga64_p_to_n_Ge64_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // d_p_to_n_p_p_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_d);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_d);
            }
#endif
            derived_rates::rate_d_p_to_n_p_p_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // He4_He4_n_to_d_Li7_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4);
            }
#endif
            derived_rates::rate_He4_He4_n_to_d_Li7_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_He4_p_to_n_B8_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_p_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_p_Be8);
            }
#endif
            derived_rates::rate_He4_He4_p_to_n_B8_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // d_He4_He4_to_p_Be9_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_d_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_d_Be8);
            }
#endif
            derived_rates::rate_d_He4_He4_to_p_Be9_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_He4_He4_to_p_B11_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_He4_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_He4_Be8);
            }
#endif
            derived_rates::rate_He4_He4_He4_to_p_B11_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_He4_p_n_to_He3_Li7_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_p_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_p_Be8);
            }
#endif
            derived_rates::rate_He4_He4_p_n_to_He3_Li7_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // He4_He4_p_n_to_p_Be9_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_p_Be8);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_p_Be8);
            }
#endif
            derived_rates::rate_He4_He4_p_n_to_p_Be9_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
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
            // Be7_to_Li7_electron_capture_reaclib
            reaclib_rates::rate_Be7_to_Li7_electron_capture_reaclib(rate_eval, tfactors);
        }

        {
            // C14_to_N14_beta_neg_reaclib
            reaclib_rates::rate_C14_to_N14_beta_neg_reaclib(rate_eval, tfactors);
        }

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
            // F17_to_O17_beta_pos_reaclib
            reaclib_rates::rate_F17_to_O17_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // F18_to_O18_beta_pos_reaclib
            reaclib_rates::rate_F18_to_O18_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Ne18_to_F18_beta_pos_reaclib
            reaclib_rates::rate_Ne18_to_F18_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Ne19_to_F19_beta_pos_reaclib
            reaclib_rates::rate_Ne19_to_F19_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Zn59_to_Cu59_beta_pos_reaclib
            reaclib_rates::rate_Zn59_to_Cu59_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Ge63_to_Ga63_beta_pos_reaclib
            reaclib_rates::rate_Ge63_to_Ga63_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // B8_to_He4_He4_beta_pos_reaclib
            reaclib_rates::rate_B8_to_He4_He4_beta_pos_reaclib(rate_eval, tfactors);
        }

        {
            // Zn59_to_p_Ni58_beta_pos_reaclib
            reaclib_rates::rate_Zn59_to_p_Ni58_beta_pos_reaclib(rate_eval, tfactors);
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
            // Na21_to_Ne21_electron_capture_weaktab
            tabular_weak_rates::rate_Na21_to_Ne21_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ne21_to_Na21_beta_neg_weaktab
            tabular_weak_rates::rate_Ne21_to_Na21_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Mg23_to_Na23_electron_capture_weaktab
            tabular_weak_rates::rate_Mg23_to_Na23_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Na23_to_Mg23_beta_neg_weaktab
            tabular_weak_rates::rate_Na23_to_Mg23_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Al25_to_Mg25_electron_capture_weaktab
            tabular_weak_rates::rate_Al25_to_Mg25_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mg25_to_Al25_beta_neg_weaktab
            tabular_weak_rates::rate_Mg25_to_Al25_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Al26_to_Mg26_electron_capture_weaktab
            tabular_weak_rates::rate_Al26_to_Mg26_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mg26_to_Al26_beta_neg_weaktab
            tabular_weak_rates::rate_Mg26_to_Al26_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // P30_to_Si30_electron_capture_weaktab
            tabular_weak_rates::rate_P30_to_Si30_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si30_to_P30_beta_neg_weaktab
            tabular_weak_rates::rate_Si30_to_P30_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // P32_to_S32_beta_neg_weaktab
            tabular_weak_rates::rate_P32_to_S32_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P32_to_Si32_electron_capture_weaktab
            tabular_weak_rates::rate_P32_to_Si32_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S32_to_P32_electron_capture_weaktab
            tabular_weak_rates::rate_S32_to_P32_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si32_to_P32_beta_neg_weaktab
            tabular_weak_rates::rate_Si32_to_P32_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Cl34_to_S34_electron_capture_weaktab
            tabular_weak_rates::rate_Cl34_to_S34_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S34_to_Cl34_beta_neg_weaktab
            tabular_weak_rates::rate_S34_to_Cl34_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Ar36_to_Cl36_electron_capture_weaktab
            tabular_weak_rates::rate_Ar36_to_Cl36_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl36_to_Ar36_beta_neg_weaktab
            tabular_weak_rates::rate_Cl36_to_Ar36_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl36_to_S36_electron_capture_weaktab
            tabular_weak_rates::rate_Cl36_to_S36_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S36_to_Cl36_beta_neg_weaktab
            tabular_weak_rates::rate_S36_to_Cl36_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Ar38_to_K38_beta_neg_weaktab
            tabular_weak_rates::rate_Ar38_to_K38_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K38_to_Ar38_electron_capture_weaktab
            tabular_weak_rates::rate_K38_to_Ar38_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar39_to_K39_beta_neg_weaktab
            tabular_weak_rates::rate_Ar39_to_K39_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K39_to_Ar39_electron_capture_weaktab
            tabular_weak_rates::rate_K39_to_Ar39_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar40_to_K40_beta_neg_weaktab
            tabular_weak_rates::rate_Ar40_to_K40_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca40_to_K40_electron_capture_weaktab
            tabular_weak_rates::rate_Ca40_to_K40_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K40_to_Ar40_electron_capture_weaktab
            tabular_weak_rates::rate_K40_to_Ar40_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K40_to_Ca40_beta_neg_weaktab
            tabular_weak_rates::rate_K40_to_Ca40_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca41_to_K41_electron_capture_weaktab
            tabular_weak_rates::rate_Ca41_to_K41_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K41_to_Ca41_beta_neg_weaktab
            tabular_weak_rates::rate_K41_to_Ca41_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca43_to_Sc43_beta_neg_weaktab
            tabular_weak_rates::rate_Ca43_to_Sc43_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc43_to_Ca43_electron_capture_weaktab
            tabular_weak_rates::rate_Sc43_to_Ca43_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca44_to_Sc44_beta_neg_weaktab
            tabular_weak_rates::rate_Ca44_to_Sc44_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc44_to_Ca44_electron_capture_weaktab
            tabular_weak_rates::rate_Sc44_to_Ca44_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc44_to_Ti44_beta_neg_weaktab
            tabular_weak_rates::rate_Sc44_to_Ti44_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti44_to_Sc44_electron_capture_weaktab
            tabular_weak_rates::rate_Ti44_to_Sc44_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co53_to_Fe53_electron_capture_weaktab
            tabular_weak_rates::rate_Co53_to_Fe53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe53_to_Co53_beta_neg_weaktab
            tabular_weak_rates::rate_Fe53_to_Co53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu57_to_Ni57_electron_capture_weaktab
            tabular_weak_rates::rate_Cu57_to_Ni57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Cu57_beta_neg_weaktab
            tabular_weak_rates::rate_Ni57_to_Cu57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca45_to_Sc45_beta_neg_weaktab
            tabular_weak_rates::rate_Ca45_to_Sc45_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc45_to_Ca45_electron_capture_weaktab
            tabular_weak_rates::rate_Sc45_to_Ca45_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc45_to_Ti45_beta_neg_weaktab
            tabular_weak_rates::rate_Sc45_to_Ti45_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti45_to_Sc45_electron_capture_weaktab
            tabular_weak_rates::rate_Ti45_to_Sc45_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca46_to_Sc46_beta_neg_weaktab
            tabular_weak_rates::rate_Ca46_to_Sc46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc46_to_Ca46_electron_capture_weaktab
            tabular_weak_rates::rate_Sc46_to_Ca46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc46_to_Ti46_beta_neg_weaktab
            tabular_weak_rates::rate_Sc46_to_Ti46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti46_to_Sc46_electron_capture_weaktab
            tabular_weak_rates::rate_Ti46_to_Sc46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti46_to_V46_beta_neg_weaktab
            tabular_weak_rates::rate_Ti46_to_V46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V46_to_Ti46_electron_capture_weaktab
            tabular_weak_rates::rate_V46_to_Ti46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca47_to_Sc47_beta_neg_weaktab
            tabular_weak_rates::rate_Ca47_to_Sc47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc47_to_Ca47_electron_capture_weaktab
            tabular_weak_rates::rate_Sc47_to_Ca47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc47_to_Ti47_beta_neg_weaktab
            tabular_weak_rates::rate_Sc47_to_Ti47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti47_to_Sc47_electron_capture_weaktab
            tabular_weak_rates::rate_Ti47_to_Sc47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti47_to_V47_beta_neg_weaktab
            tabular_weak_rates::rate_Ti47_to_V47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V47_to_Ti47_electron_capture_weaktab
            tabular_weak_rates::rate_V47_to_Ti47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca48_to_Sc48_beta_neg_weaktab
            tabular_weak_rates::rate_Ca48_to_Sc48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr48_to_V48_electron_capture_weaktab
            tabular_weak_rates::rate_Cr48_to_V48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc48_to_Ca48_electron_capture_weaktab
            tabular_weak_rates::rate_Sc48_to_Ca48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc48_to_Ti48_beta_neg_weaktab
            tabular_weak_rates::rate_Sc48_to_Ti48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti48_to_Sc48_electron_capture_weaktab
            tabular_weak_rates::rate_Ti48_to_Sc48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti48_to_V48_beta_neg_weaktab
            tabular_weak_rates::rate_Ti48_to_V48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V48_to_Cr48_beta_neg_weaktab
            tabular_weak_rates::rate_V48_to_Cr48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V48_to_Ti48_electron_capture_weaktab
            tabular_weak_rates::rate_V48_to_Ti48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr49_to_V49_electron_capture_weaktab
            tabular_weak_rates::rate_Cr49_to_V49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc49_to_Ti49_beta_neg_weaktab
            tabular_weak_rates::rate_Sc49_to_Ti49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti49_to_Sc49_electron_capture_weaktab
            tabular_weak_rates::rate_Ti49_to_Sc49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti49_to_V49_beta_neg_weaktab
            tabular_weak_rates::rate_Ti49_to_V49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V49_to_Cr49_beta_neg_weaktab
            tabular_weak_rates::rate_V49_to_Cr49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V49_to_Ti49_electron_capture_weaktab
            tabular_weak_rates::rate_V49_to_Ti49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr50_to_Mn50_beta_neg_weaktab
            tabular_weak_rates::rate_Cr50_to_Mn50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr50_to_V50_electron_capture_weaktab
            tabular_weak_rates::rate_Cr50_to_V50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn50_to_Cr50_electron_capture_weaktab
            tabular_weak_rates::rate_Mn50_to_Cr50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti50_to_V50_beta_neg_weaktab
            tabular_weak_rates::rate_Ti50_to_V50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V50_to_Cr50_beta_neg_weaktab
            tabular_weak_rates::rate_V50_to_Cr50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V50_to_Ti50_electron_capture_weaktab
            tabular_weak_rates::rate_V50_to_Ti50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr51_to_Mn51_beta_neg_weaktab
            tabular_weak_rates::rate_Cr51_to_Mn51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr51_to_V51_electron_capture_weaktab
            tabular_weak_rates::rate_Cr51_to_V51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn51_to_Cr51_electron_capture_weaktab
            tabular_weak_rates::rate_Mn51_to_Cr51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti51_to_V51_beta_neg_weaktab
            tabular_weak_rates::rate_Ti51_to_V51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V51_to_Cr51_beta_neg_weaktab
            tabular_weak_rates::rate_V51_to_Cr51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V51_to_Ti51_electron_capture_weaktab
            tabular_weak_rates::rate_V51_to_Ti51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr52_to_Mn52_beta_neg_weaktab
            tabular_weak_rates::rate_Cr52_to_Mn52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr52_to_V52_electron_capture_weaktab
            tabular_weak_rates::rate_Cr52_to_V52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe52_to_Mn52_electron_capture_weaktab
            tabular_weak_rates::rate_Fe52_to_Mn52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Cr52_electron_capture_weaktab
            tabular_weak_rates::rate_Mn52_to_Cr52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Fe52_beta_neg_weaktab
            tabular_weak_rates::rate_Mn52_to_Fe52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V52_to_Cr52_beta_neg_weaktab
            tabular_weak_rates::rate_V52_to_Cr52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr53_to_Mn53_beta_neg_weaktab
            tabular_weak_rates::rate_Cr53_to_Mn53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe53_to_Mn53_electron_capture_weaktab
            tabular_weak_rates::rate_Fe53_to_Mn53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn53_to_Cr53_electron_capture_weaktab
            tabular_weak_rates::rate_Mn53_to_Cr53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn53_to_Fe53_beta_neg_weaktab
            tabular_weak_rates::rate_Mn53_to_Fe53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co54_to_Fe54_electron_capture_weaktab
            tabular_weak_rates::rate_Co54_to_Fe54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr54_to_Mn54_beta_neg_weaktab
            tabular_weak_rates::rate_Cr54_to_Mn54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe54_to_Co54_beta_neg_weaktab
            tabular_weak_rates::rate_Fe54_to_Co54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe54_to_Mn54_electron_capture_weaktab
            tabular_weak_rates::rate_Fe54_to_Mn54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn54_to_Cr54_electron_capture_weaktab
            tabular_weak_rates::rate_Mn54_to_Cr54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn54_to_Fe54_beta_neg_weaktab
            tabular_weak_rates::rate_Mn54_to_Fe54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Fe55_to_Mn55_electron_capture_weaktab
            tabular_weak_rates::rate_Fe55_to_Mn55_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn55_to_Fe55_beta_neg_weaktab
            tabular_weak_rates::rate_Mn55_to_Fe55_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Co57_to_Fe57_electron_capture_weaktab
            tabular_weak_rates::rate_Co57_to_Fe57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co57_to_Ni57_beta_neg_weaktab
            tabular_weak_rates::rate_Co57_to_Ni57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe57_to_Co57_beta_neg_weaktab
            tabular_weak_rates::rate_Fe57_to_Co57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Co57_electron_capture_weaktab
            tabular_weak_rates::rate_Ni57_to_Co57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co58_to_Fe58_electron_capture_weaktab
            tabular_weak_rates::rate_Co58_to_Fe58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co58_to_Ni58_beta_neg_weaktab
            tabular_weak_rates::rate_Co58_to_Ni58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu58_to_Ni58_electron_capture_weaktab
            tabular_weak_rates::rate_Cu58_to_Ni58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe58_to_Co58_beta_neg_weaktab
            tabular_weak_rates::rate_Fe58_to_Co58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni58_to_Co58_electron_capture_weaktab
            tabular_weak_rates::rate_Ni58_to_Co58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni58_to_Cu58_beta_neg_weaktab
            tabular_weak_rates::rate_Ni58_to_Cu58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co59_to_Ni59_beta_neg_weaktab
            tabular_weak_rates::rate_Co59_to_Ni59_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu59_to_Ni59_electron_capture_weaktab
            tabular_weak_rates::rate_Cu59_to_Ni59_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni59_to_Co59_electron_capture_weaktab
            tabular_weak_rates::rate_Ni59_to_Co59_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni59_to_Cu59_beta_neg_weaktab
            tabular_weak_rates::rate_Ni59_to_Cu59_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu60_to_Ni60_electron_capture_weaktab
            tabular_weak_rates::rate_Cu60_to_Ni60_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu60_to_Zn60_beta_neg_weaktab
            tabular_weak_rates::rate_Cu60_to_Zn60_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni60_to_Cu60_beta_neg_weaktab
            tabular_weak_rates::rate_Ni60_to_Cu60_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn60_to_Cu60_electron_capture_weaktab
            tabular_weak_rates::rate_Zn60_to_Cu60_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu61_to_Ni61_electron_capture_weaktab
            tabular_weak_rates::rate_Cu61_to_Ni61_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu61_to_Zn61_beta_neg_weaktab
            tabular_weak_rates::rate_Cu61_to_Zn61_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni61_to_Cu61_beta_neg_weaktab
            tabular_weak_rates::rate_Ni61_to_Cu61_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn61_to_Cu61_electron_capture_weaktab
            tabular_weak_rates::rate_Zn61_to_Cu61_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu62_to_Ni62_electron_capture_weaktab
            tabular_weak_rates::rate_Cu62_to_Ni62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu62_to_Zn62_beta_neg_weaktab
            tabular_weak_rates::rate_Cu62_to_Zn62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga62_to_Zn62_electron_capture_weaktab
            tabular_weak_rates::rate_Ga62_to_Zn62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni62_to_Cu62_beta_neg_weaktab
            tabular_weak_rates::rate_Ni62_to_Cu62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn62_to_Cu62_electron_capture_weaktab
            tabular_weak_rates::rate_Zn62_to_Cu62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn62_to_Ga62_beta_neg_weaktab
            tabular_weak_rates::rate_Zn62_to_Ga62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu63_to_Ni63_electron_capture_weaktab
            tabular_weak_rates::rate_Cu63_to_Ni63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu63_to_Zn63_beta_neg_weaktab
            tabular_weak_rates::rate_Cu63_to_Zn63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga63_to_Zn63_electron_capture_weaktab
            tabular_weak_rates::rate_Ga63_to_Zn63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni63_to_Cu63_beta_neg_weaktab
            tabular_weak_rates::rate_Ni63_to_Cu63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn63_to_Cu63_electron_capture_weaktab
            tabular_weak_rates::rate_Zn63_to_Cu63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn63_to_Ga63_beta_neg_weaktab
            tabular_weak_rates::rate_Zn63_to_Ga63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu64_to_Ni64_electron_capture_weaktab
            tabular_weak_rates::rate_Cu64_to_Ni64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu64_to_Zn64_beta_neg_weaktab
            tabular_weak_rates::rate_Cu64_to_Zn64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga64_to_Ge64_beta_neg_weaktab
            tabular_weak_rates::rate_Ga64_to_Ge64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga64_to_Zn64_electron_capture_weaktab
            tabular_weak_rates::rate_Ga64_to_Zn64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ge64_to_Ga64_electron_capture_weaktab
            tabular_weak_rates::rate_Ge64_to_Ga64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni64_to_Cu64_beta_neg_weaktab
            tabular_weak_rates::rate_Ni64_to_Cu64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn64_to_Cu64_electron_capture_weaktab
            tabular_weak_rates::rate_Zn64_to_Cu64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn64_to_Ga64_beta_neg_weaktab
            tabular_weak_rates::rate_Zn64_to_Ga64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu65_to_Zn65_beta_neg_weaktab
            tabular_weak_rates::rate_Cu65_to_Zn65_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn65_to_Cu65_electron_capture_weaktab
            tabular_weak_rates::rate_Zn65_to_Cu65_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
