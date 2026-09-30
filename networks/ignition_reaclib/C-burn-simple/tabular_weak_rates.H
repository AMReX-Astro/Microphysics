#ifndef TABULAR_WEAK_RATES_H
#define TABULAR_WEAK_RATES_H

#include <cmath>

#include <AMReX.H>

#include <actual_network.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <rate_type.H>
#include <table_rates.H>

namespace tabular_weak_rates {

    using namespace Rates;
    using namespace Species;

    using namespace rate_tables;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_n_to_p_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // n --> p

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_n_to_p_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_to_p_beta_neg_weaktab) = table_values.drate_dT;
            rate_eval.denuc_weak_dY(N) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(N) * table_values.denu_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(N) * (table_values.enu + table_values.gamma);
    }


    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const amrex::Real temp,
               [[maybe_unused]] const amrex::Real rhoy,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        // Calculate tabular weak rates

        [[maybe_unused]] const amrex::Real log_temp = std::log10(temp);
        [[maybe_unused]] const amrex::Real log_rhoy = std::log10(rhoy);

        rate_eval.enuc_weak = 0.0_rt;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.denuc_weak_dT = 0.0_rt;
            for (int n = 1; n <= NumSpec; ++n) {
                rate_eval.denuc_weak_dY(n) = 0.0_rt;
            }
        }

        {
            // n_to_p_beta_neg_weaktab
            rate_n_to_p_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

    }

}
#endif
