#ifndef DERIVED_RATES_H
#define DERIVED_RATES_H

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>
#include <temperature_table_rates.H>
#include <interp_tools.H>

using namespace Rates;

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_to_He4_He4_He4_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // C12 --> 3 He4

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = He4_log_pf + He4_log_pf + He4_log_pf - C12_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dHe4_log_pf_dT9 + dHe4_log_pf_dT9 - dC12_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_He4_He4_to_C12_starlib_data::log_t9,
                                               He4_He4_He4_to_C12_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_He4_He4_to_C12_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_He4_He4_to_C12_starlib_data::log_t9,
                                                 He4_He4_He4_to_C12_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.274747609999003 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 46.74930397991297 + Q_kBT + net_log_pf + log_scor;
    _rate += 3.0 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 3.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_C12_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + C12_log_pf - He4_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dC12_log_pf_dT9 - dHe4_log_pf_dT9 - dNe20_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               C12_C12_to_He4_Ne20_starlib_data::log_t9,
                                               C12_C12_to_He4_Ne20_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_C12_C12_to_He4_Ne20_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 C12_C12_to_He4_Ne20_starlib_data::log_t9,
                                                 C12_C12_to_He4_Ne20_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -4.617016329998765 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.18785146891918028 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Na23_to_C12_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na23 + p --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = C12_log_pf + C12_log_pf - p_log_pf - Na23_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dC12_log_pf_dT9 - dp_log_pf_dT9 - dNa23_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               C12_C12_to_p_Na23_starlib_data::log_t9,
                                               C12_C12_to_p_Na23_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_C12_C12_to_p_Na23_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 C12_C12_to_p_Na23_starlib_data::log_t9,
                                                 C12_C12_to_p_Na23_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -2.2408824359984667 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += -0.039484252594465044 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_C12_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 + He4 --> C12 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + O16_log_pf - He4_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dO16_log_pf_dT9 - dHe4_log_pf_dT9 - dMg24_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               C12_O16_to_He4_Mg24_starlib_data::log_t9,
                                               C12_O16_to_He4_Mg24_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_C12_O16_to_He4_Mg24_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 C12_O16_to_He4_Mg24_starlib_data::log_t9,
                                                 C12_O16_to_He4_Mg24_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -6.771660029999435 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 1.0388581923104747 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Al27_to_C12_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Al27 + p --> C12 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + O16_log_pf - p_log_pf - Al27_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dO16_log_pf_dT9 - dp_log_pf_dT9 - dAl27_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               C12_O16_to_p_Al27_starlib_data::log_t9,
                                               C12_O16_to_p_Al27_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_C12_O16_to_p_Al27_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 C12_O16_to_p_Al27_starlib_data::log_t9,
                                                 C12_O16_to_p_Al27_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -5.170886835997453 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.43900994536895793 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Si28_to_O16_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 + He4 --> O16 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = O16_log_pf + O16_log_pf - He4_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dO16_log_pf_dT9 + dO16_log_pf_dT9 - dHe4_log_pf_dT9 - dSi28_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               O16_O16_to_He4_Si28_starlib_data::log_t9,
                                               O16_O16_to_He4_Si28_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_O16_O16_to_He4_Si28_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 O16_O16_to_He4_Si28_starlib_data::log_t9,
                                                 O16_O16_to_He4_Si28_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -9.593877029998845 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.5457157206902232 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_P31_to_O16_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // P31 + p --> O16 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real net_log_pf = O16_log_pf + O16_log_pf - p_log_pf - P31_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dO16_log_pf_dT9 + dO16_log_pf_dT9 - dp_log_pf_dT9 - dP31_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               O16_O16_to_p_P31_starlib_data::log_t9,
                                               O16_O16_to_p_P31_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_O16_O16_to_p_P31_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 O16_O16_to_p_P31_starlib_data::log_t9,
                                                 O16_O16_to_p_P31_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.677569135998522 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 1.0683848103569353 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_O16_to_He4_N13_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 + p --> He4 + N13

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real N13_log_pf, dN13_log_pf_dT9;

    // setting N13 log(partition function) to 0.0 by default, independent of T
    N13_log_pf = 0.0_rt;
    dN13_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = He4_log_pf + N13_log_pf - p_log_pf - O16_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dN13_log_pf_dT9 - dp_log_pf_dT9 - dO16_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_N13_to_p_O16_starlib_data::log_t9,
                                               He4_N13_to_p_O16_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_N13_to_p_O16_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_N13_to_p_O16_starlib_data::log_t9,
                                                 He4_N13_to_p_O16_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -5.218426906001445 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 1.7520206134245597 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_N13_to_p_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N13 --> p + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real N13_log_pf, dN13_log_pf_dT9;

    // setting N13 log(partition function) to 0.0 by default, independent of T
    N13_log_pf = 0.0_rt;
    dN13_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + C12_log_pf - N13_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dC12_log_pf_dT9 - dN13_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_C12_to_N13_starlib_data::log_t9,
                                               p_C12_to_N13_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_C12_to_N13_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_C12_to_N13_starlib_data::log_t9,
                                                 p_C12_to_N13_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -1.9434910639993015 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 22.910928683693466 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_O16_to_He4_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 --> He4 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + C12_log_pf - O16_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dC12_log_pf_dT9 - dO16_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_C12_to_O16_starlib_data::log_t9,
                                               He4_C12_to_O16_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_C12_to_O16_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_C12_to_O16_starlib_data::log_t9,
                                                 He4_C12_to_O16_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.161917969999712 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.662949297118026 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ne20_to_He4_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 --> He4 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + O16_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dO16_log_pf_dT9 - dNe20_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_O16_to_Ne20_starlib_data::log_t9,
                                               He4_O16_to_Ne20_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_O16_to_Ne20_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_O16_to_Ne20_starlib_data::log_t9,
                                                 He4_O16_to_Ne20_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -4.729845970000538 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.759112801502127 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_p_Na23_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Na23_log_pf - He4_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNa23_log_pf_dT9 - dHe4_log_pf_dT9 - dNe20_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Na23_to_He4_Ne20_starlib_data::log_t9,
                                               p_Na23_to_He4_Ne20_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Na23_to_He4_Ne20_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Na23_to_He4_Ne20_starlib_data::log_t9,
                                                 p_Na23_to_He4_Ne20_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -2.3761338939988508 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.22733572151364528 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mg24_to_p_Na23_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Na23_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNa23_log_pf_dT9 - dMg24_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Na23_to_Mg24_starlib_data::log_t9,
                                               p_Na23_to_Mg24_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Na23_to_Mg24_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Na23_to_Mg24_starlib_data::log_t9,
                                                 p_Na23_to_Mg24_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -11.692695564000061 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.04814456146302 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mg24_to_He4_Ne20_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ne20_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNe20_log_pf_dT9 - dMg24_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Ne20_to_Mg24_starlib_data::log_t9,
                                               He4_Ne20_to_Mg24_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Ne20_to_Mg24_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Ne20_to_Mg24_starlib_data::log_t9,
                                                 He4_Ne20_to_Mg24_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -9.31656167000121 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.820808839949375 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_p_Al27_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 + He4 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + Al27_log_pf - He4_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dAl27_log_pf_dT9 - dHe4_log_pf_dT9 - dMg24_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Al27_to_He4_Mg24_starlib_data::log_t9,
                                               p_Al27_to_He4_Mg24_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Al27_to_He4_Mg24_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Al27_to_He4_Mg24_starlib_data::log_t9,
                                                 p_Al27_to_He4_Mg24_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -1.6007731940019818 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.5998482469415165 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Si28_to_p_Al27_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Al27_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dAl27_log_pf_dT9 - dSi28_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Al27_to_Si28_starlib_data::log_t9,
                                               p_Al27_to_Si28_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Al27_to_Si28_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Al27_to_Si28_starlib_data::log_t9,
                                                 p_Al27_to_Si28_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -11.584908164001932 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.462802252999236 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Si28_to_He4_Mg24_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> He4 + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mg24_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMg24_log_pf_dT9 - dSi28_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Mg24_to_Si28_starlib_data::log_t9,
                                               He4_Mg24_to_Si28_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Mg24_to_Si28_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Mg24_to_Si28_starlib_data::log_t9,
                                                 He4_Mg24_to_Si28_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -9.98413496999995 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.86295400605772 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Si28_to_p_P31_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 + He4 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + P31_log_pf - He4_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dP31_log_pf_dT9 - dHe4_log_pf_dT9 - dSi28_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_P31_to_He4_Si28_starlib_data::log_t9,
                                               p_P31_to_He4_Si28_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_P31_to_He4_Si28_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_P31_to_He4_Si28_starlib_data::log_t9,
                                                 p_P31_to_He4_Si28_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -1.9163078940003233 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += -0.5226690896667122 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_S32_to_p_P31_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + P31_log_pf - S32_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dP31_log_pf_dT9 - dS32_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_P31_to_S32_starlib_data::log_t9,
                                               p_P31_to_S32_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_P31_to_S32_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_P31_to_S32_starlib_data::log_t9,
                                                 p_P31_to_S32_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.863963764002593 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.370833422303182 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_S32_to_He4_Si28_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Si28_log_pf - S32_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSi28_log_pf_dT9 - dS32_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Si28_to_S32_starlib_data::log_t9,
                                               He4_Si28_to_S32_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Si28_to_S32_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Si28_to_S32_starlib_data::log_t9,
                                                 He4_Si28_to_S32_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -6.9476558700022695 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.893502511969896 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mn51_to_He4_Cr48_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn51 + p --> He4 + Cr48

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cr48_log_pf - p_log_pf - Mn51_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr48_log_pf_dT9 - dp_log_pf_dT9 - dMn51_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Cr48_to_p_Mn51_starlib_data::log_t9,
                                               He4_Cr48_to_p_Mn51_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Cr48_to_p_Mn51_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Cr48_to_p_Mn51_starlib_data::log_t9,
                                                 He4_Cr48_to_p_Mn51_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -0.5581448060026141 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += -0.5141173297653887 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co55_to_He4_Fe52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 + p --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe52_log_pf - p_log_pf - Co55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe52_log_pf_dT9 - dp_log_pf_dT9 - dCo55_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Fe52_to_p_Co55_starlib_data::log_t9,
                                               He4_Fe52_to_p_Co55_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Fe52_to_p_Co55_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Fe52_to_p_Co55_starlib_data::log_t9,
                                                 He4_Fe52_to_p_Co55_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -0.8338448059985103 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += -0.7949991911878357 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe52_to_p_Mn51_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe52 --> p + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Mn51_log_pf - Fe52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMn51_log_pf_dT9 - dFe52_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Mn51_to_Fe52_starlib_data::log_t9,
                                               p_Mn51_to_Fe52_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Mn51_to_Fe52_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Mn51_to_Fe52_starlib_data::log_t9,
                                                 p_Mn51_to_Fe52_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.377871064004735 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.487831322280396 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe52_to_He4_Cr48_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe52 --> He4 + Cr48

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cr48_log_pf - Fe52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr48_log_pf_dT9 - dFe52_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Cr48_to_Fe52_starlib_data::log_t9,
                                               He4_Cr48_to_Fe52_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Cr48_to_Fe52_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Cr48_to_Fe52_starlib_data::log_t9,
                                                 He4_Cr48_to_Fe52_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.9360158700054875 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.973713992515005 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co55_to_He4_Mn51_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 --> He4 + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mn51_log_pf - Co55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn51_log_pf_dT9 - dCo55_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Mn51_to_Co55_starlib_data::log_t9,
                                               He4_Mn51_to_Co55_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Mn51_to_Co55_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Mn51_to_Co55_starlib_data::log_t9,
                                                 He4_Mn51_to_Co55_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.211715870003866 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.692832131092562 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni56_to_p_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> p + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Co55_log_pf - Ni56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCo55_log_pf_dT9 - dNi56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Co55_to_Ni56_starlib_data::log_t9,
                                               p_Co55_to_Ni56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Co55_to_Ni56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Co55_to_Ni56_starlib_data::log_t9,
                                                 p_Co55_to_Ni56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.166571063993911 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.777604041925606 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni56_to_He4_Fe52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe52_log_pf - Ni56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe52_log_pf_dT9 - dNi56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Fe52_to_Ni56_starlib_data::log_t9,
                                               He4_Fe52_to_Ni56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Fe52_to_Ni56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Fe52_to_Ni56_starlib_data::log_t9,
                                                 He4_Fe52_to_Ni56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.000415869992421 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.98260485073777 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ni56_to_p_Cu59_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 + He4 --> p + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + Cu59_log_pf - He4_log_pf - Ni56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCu59_log_pf_dT9 - dHe4_log_pf_dT9 - dNi56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Cu59_to_He4_Ni56_starlib_data::log_t9,
                                               p_Cu59_to_He4_Ni56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Cu59_to_He4_Ni56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Cu59_to_He4_Ni56_starlib_data::log_t9,
                                                 p_Cu59_to_He4_Ni56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -2.4131551940024294 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.09615846250448774 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cu59_to_He4_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cu59 --> He4 + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Co55_log_pf - Cu59_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCo55_log_pf_dT9 - dCu59_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Co55_to_Cu59_starlib_data::log_t9,
                                               He4_Co55_to_Cu59_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Co55_to_Cu59_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Co55_to_Cu59_starlib_data::log_t9,
                                                 He4_Co55_to_Cu59_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -4.7534158699933435 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.68144557942112 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Zn60_to_p_Cu59_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Zn60 --> p + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Zn60_log_pf, dZn60_log_pf_dT9;

    // interpolating Zn60 partition function
    get_partition_function_cached(Zn60, tfactors.T9, pf_cache, Zn60_log_pf, dZn60_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Cu59_log_pf - Zn60_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCu59_log_pf_dT9 - dZn60_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Cu59_to_Zn60_starlib_data::log_t9,
                                               p_Cu59_to_Zn60_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Cu59_to_Zn60_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Cu59_to_Zn60_starlib_data::log_t9,
                                                 p_Cu59_to_Zn60_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -5.104971064000085 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.08618638341413 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Zn60_to_He4_Ni56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Zn60 --> He4 + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Zn60_log_pf, dZn60_log_pf_dT9;

    // interpolating Zn60 partition function
    get_partition_function_cached(Zn60, tfactors.T9, pf_cache, Zn60_log_pf, dZn60_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ni56_log_pf - Zn60_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNi56_log_pf_dT9 - dZn60_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Ni56_to_Zn60_starlib_data::log_t9,
                                               He4_Ni56_to_Zn60_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Ni56_to_Zn60_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Ni56_to_Zn60_starlib_data::log_t9,
                                                 He4_Ni56_to_Zn60_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -2.6918158700001382 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.990027920909643 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mg24_to_C12_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + C12_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dC12_log_pf_dT9 - dMg24_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               C12_C12_to_n_Mg23_starlib_data::log_t9,
                                               C12_C12_to_n_Mg23_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_C12_C12_to_n_Mg23_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 C12_C12_to_n_Mg23_starlib_data::log_t9,
                                                 C12_C12_to_n_Mg23_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -13.933577999999976 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.008660308868553 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_S32_to_O16_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> O16 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = O16_log_pf + O16_log_pf - S32_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dO16_log_pf_dT9 + dO16_log_pf_dT9 - dS32_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               O16_O16_to_n_S31_starlib_data::log_t9,
                                               O16_O16_to_n_S31_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_O16_O16_to_n_S31_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 O16_O16_to_n_S31_starlib_data::log_t9,
                                                 O16_O16_to_n_S31_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -16.541532900001116 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.439218232660117 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Si28_to_C12_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> C12 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + O16_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dO16_log_pf_dT9 - dSi28_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               C12_O16_to_n_Si27_starlib_data::log_t9,
                                               C12_O16_to_n_Si27_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_C12_O16_to_n_Si27_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 C12_O16_to_n_Si27_starlib_data::log_t9,
                                                 C12_O16_to_n_Si27_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -16.755794999999384 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.901812198368194 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe54_to_He4_Mn51_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 + p --> He4 + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mn51_log_pf - p_log_pf - Fe54_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn51_log_pf_dT9 - dp_log_pf_dT9 - dFe54_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Mn51_to_p_Fe54_starlib_data::log_t9,
                                               He4_Mn51_to_p_Fe54_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Mn51_to_p_Fe54_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Mn51_to_p_Fe54_starlib_data::log_t9,
                                                 He4_Mn51_to_p_Fe54_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -3.1473448060034297 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 3.074713682988042 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe56_to_n_Co56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe56 + p --> n + Co56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co56_log_pf - p_log_pf - Fe56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo56_log_pf_dT9 - dp_log_pf_dT9 - dFe56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Co56_to_p_Fe56_starlib_data::log_t9,
                                               n_Co56_to_p_Fe56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Co56_to_p_Fe56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Co56_to_p_Fe56_starlib_data::log_t9,
                                                 n_Co56_to_p_Fe56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -5.349007036003482 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 2.1995262241642837 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_p_Co57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 + He4 --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Co57_log_pf - He4_log_pf - Fe54_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCo57_log_pf_dT9 - dHe4_log_pf_dT9 - dFe54_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Co57_to_He4_Fe54_starlib_data::log_t9,
                                               p_Co57_to_He4_Fe54_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Co57_to_He4_Fe54_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Co57_to_He4_Fe54_starlib_data::log_t9,
                                                 p_Co57_to_He4_Fe54_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -1.772955193998673 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.7921095556653313 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ni58_to_He4_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 + p --> He4 + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Co55_log_pf - p_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCo55_log_pf_dT9 - dp_log_pf_dT9 - dNi58_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Co55_to_p_Ni58_starlib_data::log_t9,
                                               He4_Co55_to_p_Ni58_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Co55_to_p_Ni58_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Co55_to_p_Ni58_starlib_data::log_t9,
                                                 He4_Co55_to_p_Ni58_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -1.3348448059966436 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 3.3683618674820073 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co55_to_p_Fe54_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 --> p + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Fe54_log_pf - Co55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dFe54_log_pf_dT9 - dCo55_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Fe54_to_Co55_starlib_data::log_t9,
                                               p_Fe54_to_Co55_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Fe54_to_Co55_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Fe54_to_Co55_starlib_data::log_t9,
                                                 p_Fe54_to_Co55_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -5.0643710640004365 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 21.618118448104518 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co56_to_n_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 --> n + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co55_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo55_log_pf_dT9 - dCo56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Co55_to_Co56_starlib_data::log_t9,
                                               n_Co55_to_Co56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Co55_to_Co56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Co55_to_Co56_starlib_data::log_t9,
                                                 n_Co55_to_Co56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -10.081818100002103 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 23.582564333711023 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co57_to_n_Co56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co57 --> n + Co56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co56_log_pf - Co57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo56_log_pf_dT9 - dCo57_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Co56_to_Co57_starlib_data::log_t9,
                                               n_Co56_to_Co57_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Co56_to_Co57_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Co56_to_Co57_starlib_data::log_t9,
                                                 n_Co56_to_Co57_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -11.376518100002114 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 23.81866809654295 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co57_to_p_Fe56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co57 --> p + Fe56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Fe56_log_pf - Co57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dFe56_log_pf_dT9 - dCo57_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Fe56_to_Co57_starlib_data::log_t9,
                                               p_Fe56_to_Co57_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Fe56_to_Co57_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Fe56_to_Co57_starlib_data::log_t9,
                                                 p_Fe56_to_Co57_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -6.027511063998631 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 21.619141872378666 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co56_to_n_Ni56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 + p --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni56_log_pf - p_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi56_log_pf_dT9 - dp_log_pf_dT9 - dCo56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Ni56_to_p_Co56_starlib_data::log_t9,
                                               n_Ni56_to_p_Co56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Ni56_to_p_Co56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Ni56_to_p_Co56_starlib_data::log_t9,
                                                 n_Ni56_to_p_Co56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -2.915247036008191 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += -2.195039708214582 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni58_to_p_Co57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Co57_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCo57_log_pf_dT9 - dNi58_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Co57_to_Ni58_starlib_data::log_t9,
                                               p_Co57_to_Ni58_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Co57_to_Ni58_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Co57_to_Ni58_starlib_data::log_t9,
                                                 p_Co57_to_Ni58_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.172171063995753 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.778589871251857 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni58_to_He4_Fe54_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> He4 + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe54_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe54_log_pf_dT9 - dNi58_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Fe54_to_Ni58_starlib_data::log_t9,
                                               He4_Fe54_to_Ni58_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Fe54_to_Ni58_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Fe54_to_Ni58_starlib_data::log_t9,
                                                 He4_Fe54_to_Ni58_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -6.39921586999708 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.986480315586526 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Co56_to_n_Cu59_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 + He4 --> n + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Cu59_log_pf - He4_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCu59_log_pf_dT9 - dHe4_log_pf_dT9 - dCo56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Cu59_to_He4_Co56_starlib_data::log_t9,
                                               n_Cu59_to_He4_Co56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Cu59_to_He4_Co56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Cu59_to_He4_Co56_starlib_data::log_t9,
                                                 n_Cu59_to_He4_Co56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -5.328402230008759 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += -2.098881245710094 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cu59_to_p_Ni58_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cu59 --> p + Ni58

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ni58_log_pf - Cu59_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNi58_log_pf_dT9 - dCu59_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Ni58_to_Cu59_starlib_data::log_t9,
                                               p_Ni58_to_Cu59_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Ni58_to_Cu59_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Ni58_to_Cu59_starlib_data::log_t9,
                                                 p_Ni58_to_Cu59_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -3.4185710639967 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 22.31308371193911 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_S32_to_p_Cl35_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 + He4 --> p + Cl35

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + Cl35_log_pf - He4_log_pf - S32_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCl35_log_pf_dT9 - dHe4_log_pf_dT9 - dS32_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Cl35_to_He4_S32_starlib_data::log_t9,
                                               p_Cl35_to_He4_S32_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Cl35_to_He4_S32_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Cl35_to_He4_S32_starlib_data::log_t9,
                                                 p_Cl35_to_He4_S32_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -1.8660622939993667 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.15223526553037453 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ar36_to_He4_S32_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 --> He4 + S32

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + S32_log_pf - Ar36_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dS32_log_pf_dT9 - dAr36_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_S32_to_Ar36_starlib_data::log_t9,
                                               He4_S32_to_Ar36_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_S32_to_Ar36_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_S32_to_Ar36_starlib_data::log_t9,
                                                 He4_S32_to_Ar36_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -6.640920770003172 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.917075071591125 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ar36_to_p_Cl35_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 --> p + Cl35

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Cl35_log_pf - Ar36_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCl35_log_pf_dT9 - dAr36_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Cl35_to_Ar36_starlib_data::log_t9,
                                               p_Cl35_to_Ar36_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Cl35_to_Ar36_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Cl35_to_Ar36_starlib_data::log_t9,
                                                 p_Cl35_to_Ar36_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.506983063999229 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.0693103371215 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ar36_to_p_K39_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 + He4 --> p + K39

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real K39_log_pf, dK39_log_pf_dT9;

    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + K39_log_pf - He4_log_pf - Ar36_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dK39_log_pf_dT9 - dHe4_log_pf_dT9 - dAr36_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_K39_to_He4_Ar36_starlib_data::log_t9,
                                               p_K39_to_He4_Ar36_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_K39_to_He4_Ar36_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_K39_to_He4_Ar36_starlib_data::log_t9,
                                                 p_K39_to_He4_Ar36_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -1.288402194006114 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.13785123272943048 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ca40_to_He4_Ar36_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 --> He4 + Ar36

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ar36_log_pf - Ca40_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAr36_log_pf_dT9 - dCa40_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Ar36_to_Ca40_starlib_data::log_t9,
                                               He4_Ar36_to_Ca40_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Ar36_to_Ca40_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Ar36_to_Ca40_starlib_data::log_t9,
                                                 He4_Ar36_to_Ca40_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.039775869997547 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.935721252509254 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ca40_to_p_K39_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 --> p + K39

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real K39_log_pf, dK39_log_pf_dT9;

    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + K39_log_pf - Ca40_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dK39_log_pf_dT9 - dCa40_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_K39_to_Ca40_starlib_data::log_t9,
                                               p_K39_to_Ca40_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_K39_to_Ca40_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_K39_to_Ca40_starlib_data::log_t9,
                                                 p_K39_to_Ca40_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.32817806400366 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.07357248523868 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ca40_to_p_Sc43_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 + He4 --> p + Sc43

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Sc43_log_pf - He4_log_pf - Ca40_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSc43_log_pf_dT9 - dHe4_log_pf_dT9 - dCa40_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Sc43_to_He4_Ca40_starlib_data::log_t9,
                                               p_Sc43_to_He4_Ca40_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Sc43_to_He4_Ca40_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Sc43_to_He4_Ca40_starlib_data::log_t9,
                                                 p_Sc43_to_He4_Ca40_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -3.522357194003101 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.819567557884143 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ti44_to_He4_Ca40_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti44 --> He4 + Ca40

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ca40_log_pf - Ti44_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCa40_log_pf_dT9 - dTi44_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Ca40_to_Ti44_starlib_data::log_t9,
                                               He4_Ca40_to_Ti44_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Ca40_to_Ti44_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Ca40_to_Ti44_starlib_data::log_t9,
                                                 He4_Ca40_to_Ti44_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -5.127113869992141 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.950664648348603 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ti44_to_p_Sc43_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti44 --> p + Sc43

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Sc43_log_pf - Ti44_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSc43_log_pf_dT9 - dTi44_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_Sc43_to_Ti44_starlib_data::log_t9,
                                               p_Sc43_to_Ti44_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_Sc43_to_Ti44_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_Sc43_to_Ti44_starlib_data::log_t9,
                                                 p_Sc43_to_Ti44_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.649471063995243 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.77023220623275 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ti44_to_p_V47_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti44 + He4 --> p + V47

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real V47_log_pf, dV47_log_pf_dT9;

    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + V47_log_pf - He4_log_pf - Ti44_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dV47_log_pf_dT9 - dHe4_log_pf_dT9 - dTi44_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_V47_to_He4_Ti44_starlib_data::log_t9,
                                               p_V47_to_He4_Ti44_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_V47_to_He4_Ti44_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_V47_to_He4_Ti44_starlib_data::log_t9,
                                                 p_V47_to_He4_Ti44_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -0.4055851939957223 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 0.1166883815247219 + Q_kBT + net_log_pf + log_scor;
    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 0.0 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cr48_to_He4_Ti44_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr48 --> He4 + Ti44

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ti44_log_pf - Cr48_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dTi44_log_pf_dT9 - dCr48_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               He4_Ti44_to_Cr48_starlib_data::log_t9,
                                               He4_Ti44_to_Cr48_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_He4_Ti44_to_Cr48_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 He4_Ti44_to_Cr48_starlib_data::log_t9,
                                                 He4_Ti44_to_Cr48_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -7.6973158699968405 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 24.963264658324956 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cr48_to_p_V47_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr48 --> p + V47

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real V47_log_pf, dV47_log_pf_dT9;

    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + V47_log_pf - Cr48_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dV47_log_pf_dT9 - dCr48_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               p_V47_to_Cr48_starlib_data::log_t9,
                                               p_V47_to_Cr48_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_p_V47_to_Cr48_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 p_V47_to_Cr48_starlib_data::log_t9,
                                                 p_V47_to_Cr48_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -8.102901063992563 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.079953039849677 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe54_to_n_Fe53_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 --> n + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe53_log_pf - Fe54_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe53_log_pf_dT9 - dFe54_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Fe53_to_Fe54_starlib_data::log_t9,
                                               n_Fe53_to_Fe54_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Fe53_to_Fe54_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Fe53_to_Fe54_starlib_data::log_t9,
                                                 n_Fe53_to_Fe54_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -13.378418100002976 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.778945367634943 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe53_to_n_Fe52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe53 --> n + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe52_log_pf - Fe53_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe52_log_pf_dT9 - dFe53_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Fe52_to_Fe53_starlib_data::log_t9,
                                               n_Fe52_to_Fe53_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Fe52_to_Fe53_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Fe52_to_Fe53_starlib_data::log_t9,
                                                 n_Fe52_to_Fe53_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -10.686718099998421 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 21.619397427096526 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe56_to_n_Fe55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe56 --> n + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe55_log_pf - Fe56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe55_log_pf_dT9 - dFe56_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Fe55_to_Fe56_starlib_data::log_t9,
                                               n_Fe55_to_Fe56_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Fe55_to_Fe56_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Fe55_to_Fe56_starlib_data::log_t9,
                                                 n_Fe55_to_Fe56_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -11.197078100005754 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.086692234523206 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe55_to_n_Fe54_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe55 --> n + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe54_log_pf - Fe55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe54_log_pf_dT9 - dFe55_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Fe54_to_Fe55_starlib_data::log_t9,
                                               n_Fe54_to_Fe55_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Fe54_to_Fe55_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Fe54_to_Fe55_starlib_data::log_t9,
                                                 n_Fe54_to_Fe55_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -9.298118100000266 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 22.31351677145662 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni58_to_n_Ni57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> n + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni57_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi57_log_pf_dT9 - dNi58_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Ni57_to_Ni58_starlib_data::log_t9,
                                               n_Ni57_to_Ni58_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Ni57_to_Ni58_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Ni57_to_Ni58_starlib_data::log_t9,
                                                 n_Ni57_to_Ni58_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -12.216218100000214 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 25.087678975583227 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni57_to_n_Ni56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni57 --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni56_log_pf - Ni57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi56_log_pf_dT9 - dNi57_log_pf_dT9;

    auto [_rate, _drate_dT] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                               tfactors.lnT9,
                                               n_Ni56_to_Ni57_starlib_data::log_t9,
                                               n_Ni56_to_Ni57_starlib_data::log_rate);

    auto p = Rates::get_p_random<k_n_Ni56_to_Ni57_starlib>();
    auto [_sigma, _dsigma_dlogT9] = interp_net::monotone_1d_interp<do_T_derivatives>(
                                                 tfactors.lnT9,
                                                 n_Ni56_to_Ni57_starlib_data::log_t9,
                                                 n_Ni56_to_Ni57_starlib_data::sigma_rate);

    // Add rate uncertainty
    _rate += p * _sigma;

    // Apply Equilibrium Ratio
    constexpr amrex::Real Q_kBGK = -10.247718100005844 * 1.0e-9_rt / C::k_MeV;
    amrex::Real Q_kBT = Q_kBGK * tfactors.T9i;
    _rate += 22.314539283996996 + Q_kBT + net_log_pf + log_scor;
    _rate += 1.5 * tfactors.lnT9;

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    _rate = std::max(_rate, -230.0);
    rate = std::exp(_rate);

    // we found dlog(rate)/dlog(T9)
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        // Add rate uncertainty term
        _drate_dT += p * _dsigma_dlogT9;

        // Convert to dlog(rate)/dT9 first
        _drate_dT = (_drate_dT + 1.5 - Q_kBT) * tfactors.T9i + net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
        drate_dT = rate * _drate_dT * 1.0e-9_rt;
    }

}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_derived_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::ignore_unused(tfactors);
    amrex::ignore_unused(rate_eval);

    [[maybe_unused]] amrex::Real rate;
    [[maybe_unused]] amrex::Real drate_dT;

    part_fun::pf_cache_t pf_cache{};

    pf_cache.index_temp_array_1 = interp_net::find_index(tfactors.T9, part_fun::temp_array_1);

    {
        // C12_to_He4_He4_He4_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_C12_to_He4_He4_He4_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_C12_to_He4_He4_He4_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4_derived) = drate_dT;
        }
    }

    {
        // He4_Ne20_to_C12_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ne20);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
        }
#endif
        rate_He4_Ne20_to_C12_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = drate_dT;
        }
    }

    {
        // p_Na23_to_C12_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Na23);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
        }
#endif
        rate_p_Na23_to_C12_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Na23_to_C12_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12_derived) = drate_dT;
        }
    }

    {
        // He4_Mg24_to_C12_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mg24);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
        }
#endif
        rate_He4_Mg24_to_C12_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Mg24_to_C12_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16_derived) = drate_dT;
        }
    }

    {
        // p_Al27_to_C12_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Al27);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
        }
#endif
        rate_p_Al27_to_C12_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Al27_to_C12_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Al27_to_C12_O16_derived) = drate_dT;
        }
    }

    {
        // He4_Si28_to_O16_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Si28);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
        }
#endif
        rate_He4_Si28_to_O16_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Si28_to_O16_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Si28_to_O16_O16_derived) = drate_dT;
        }
    }

    {
        // p_P31_to_O16_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_P31);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
        }
#endif
        rate_p_P31_to_O16_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_P31_to_O16_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_P31_to_O16_O16_derived) = drate_dT;
        }
    }

    {
        // p_O16_to_He4_N13_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O16);
        }
#endif
        rate_p_O16_to_He4_N13_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = drate_dT;
        }
    }

    {
        // N13_to_p_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_N13_to_p_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_N13_to_p_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_N13_to_p_C12_derived) = drate_dT;
        }
    }

    {
        // O16_to_He4_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_O16_to_He4_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;
        }
    }

    {
        // Ne20_to_He4_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ne20_to_He4_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ne20_to_He4_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_derived) = drate_dT;
        }
    }

    {
        // He4_Ne20_to_p_Na23_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ne20);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
        }
#endif
        rate_He4_Ne20_to_p_Na23_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ne20_to_p_Na23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23_derived) = drate_dT;
        }
    }

    {
        // Mg24_to_p_Na23_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mg24_to_p_Na23_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mg24_to_p_Na23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived) = drate_dT;
        }
    }

    {
        // Mg24_to_He4_Ne20_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mg24_to_He4_Ne20_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;
        }
    }

    {
        // He4_Mg24_to_p_Al27_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mg24);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
        }
#endif
        rate_He4_Mg24_to_p_Al27_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_derived) = drate_dT;
        }
    }

    {
        // Si28_to_p_Al27_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Si28_to_p_Al27_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Si28_to_p_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_derived) = drate_dT;
        }
    }

    {
        // Si28_to_He4_Mg24_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Si28_to_He4_Mg24_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Si28_to_He4_Mg24_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_derived) = drate_dT;
        }
    }

    {
        // He4_Si28_to_p_P31_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Si28);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
        }
#endif
        rate_He4_Si28_to_p_P31_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Si28_to_p_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_derived) = drate_dT;
        }
    }

    {
        // S32_to_p_P31_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_S32_to_p_P31_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_S32_to_p_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_p_P31_derived) = drate_dT;
        }
    }

    {
        // S32_to_He4_Si28_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_S32_to_He4_Si28_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_S32_to_He4_Si28_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_derived) = drate_dT;
        }
    }

    {
        // p_Mn51_to_He4_Cr48_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mn51);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
        }
#endif
        rate_p_Mn51_to_He4_Cr48_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_derived) = drate_dT;
        }
    }

    {
        // p_Co55_to_He4_Fe52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
        }
#endif
        rate_p_Co55_to_He4_Fe52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived) = drate_dT;
        }
    }

    {
        // Fe52_to_p_Mn51_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe52_to_p_Mn51_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe52_to_p_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_derived) = drate_dT;
        }
    }

    {
        // Fe52_to_He4_Cr48_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe52_to_He4_Cr48_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe52_to_He4_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_derived) = drate_dT;
        }
    }

    {
        // Co55_to_He4_Mn51_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co55_to_He4_Mn51_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co55_to_He4_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51_derived) = drate_dT;
        }
    }

    {
        // Ni56_to_p_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni56_to_p_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni56_to_p_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_derived) = drate_dT;
        }
    }

    {
        // Ni56_to_He4_Fe52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni56_to_He4_Fe52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni56_to_He4_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_derived) = drate_dT;
        }
    }

    {
        // He4_Ni56_to_p_Cu59_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ni56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
        }
#endif
        rate_He4_Ni56_to_p_Cu59_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ni56_to_p_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ni56_to_p_Cu59_derived) = drate_dT;
        }
    }

    {
        // Cu59_to_He4_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cu59_to_He4_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cu59_to_He4_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_He4_Co55_derived) = drate_dT;
        }
    }

    {
        // Zn60_to_p_Cu59_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Zn60_to_p_Cu59_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Zn60_to_p_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_p_Cu59_derived) = drate_dT;
        }
    }

    {
        // Zn60_to_He4_Ni56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Zn60_to_He4_Ni56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Zn60_to_He4_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_He4_Ni56_derived) = drate_dT;
        }
    }

    {
        // Mg24_to_C12_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mg24_to_C12_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mg24_to_C12_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_C12_C12_derived) = drate_dT;
        }
    }

    {
        // S32_to_O16_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_S32_to_O16_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_S32_to_O16_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_O16_O16_derived) = drate_dT;
        }
    }

    {
        // Si28_to_C12_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Si28_to_C12_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Si28_to_C12_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_C12_O16_derived) = drate_dT;
        }
    }

    {
        // p_Fe54_to_He4_Mn51_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe54);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
        }
#endif
        rate_p_Fe54_to_He4_Mn51_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Fe54_to_He4_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe54_to_He4_Mn51_derived) = drate_dT;
        }
    }

    {
        // p_Fe56_to_n_Co56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
        }
#endif
        rate_p_Fe56_to_n_Co56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Fe56_to_n_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe56_to_n_Co56_derived) = drate_dT;
        }
    }

    {
        // He4_Fe54_to_p_Co57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe54);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
        }
#endif
        rate_He4_Fe54_to_p_Co57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Fe54_to_p_Co57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_Co57_derived) = drate_dT;
        }
    }

    {
        // p_Ni58_to_He4_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ni58);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
        }
#endif
        rate_p_Ni58_to_He4_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Ni58_to_He4_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ni58_to_He4_Co55_derived) = drate_dT;
        }
    }

    {
        // Co55_to_p_Fe54_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co55_to_p_Fe54_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co55_to_p_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_p_Fe54_derived) = drate_dT;
        }
    }

    {
        // Co56_to_n_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co56_to_n_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co56_to_n_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_n_Co55_derived) = drate_dT;
        }
    }

    {
        // Co57_to_n_Co56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co57_to_n_Co56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co57_to_n_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_n_Co56_derived) = drate_dT;
        }
    }

    {
        // Co57_to_p_Fe56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co57_to_p_Fe56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co57_to_p_Fe56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_p_Fe56_derived) = drate_dT;
        }
    }

    {
        // p_Co56_to_n_Ni56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
        }
#endif
        rate_p_Co56_to_n_Ni56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Co56_to_n_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co56_to_n_Ni56_derived) = drate_dT;
        }
    }

    {
        // Ni58_to_p_Co57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni58_to_p_Co57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni58_to_p_Co57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_p_Co57_derived) = drate_dT;
        }
    }

    {
        // Ni58_to_He4_Fe54_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni58_to_He4_Fe54_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni58_to_He4_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_He4_Fe54_derived) = drate_dT;
        }
    }

    {
        // He4_Co56_to_n_Cu59_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Co56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co56);
        }
#endif
        rate_He4_Co56_to_n_Cu59_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Co56_to_n_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Co56_to_n_Cu59_derived) = drate_dT;
        }
    }

    {
        // Cu59_to_p_Ni58_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cu59_to_p_Ni58_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cu59_to_p_Ni58_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_p_Ni58_derived) = drate_dT;
        }
    }

    {
        // He4_S32_to_p_Cl35_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_S32);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
        }
#endif
        rate_He4_S32_to_p_Cl35_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_S32_to_p_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_derived) = drate_dT;
        }
    }

    {
        // Ar36_to_He4_S32_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ar36_to_He4_S32_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ar36_to_He4_S32_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_derived) = drate_dT;
        }
    }

    {
        // Ar36_to_p_Cl35_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ar36_to_p_Cl35_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ar36_to_p_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_derived) = drate_dT;
        }
    }

    {
        // He4_Ar36_to_p_K39_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ar36);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
        }
#endif
        rate_He4_Ar36_to_p_K39_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ar36_to_p_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_derived) = drate_dT;
        }
    }

    {
        // Ca40_to_He4_Ar36_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ca40_to_He4_Ar36_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ca40_to_He4_Ar36_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_derived) = drate_dT;
        }
    }

    {
        // Ca40_to_p_K39_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ca40_to_p_K39_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ca40_to_p_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_derived) = drate_dT;
        }
    }

    {
        // He4_Ca40_to_p_Sc43_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ca40);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
        }
#endif
        rate_He4_Ca40_to_p_Sc43_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_derived) = drate_dT;
        }
    }

    {
        // Ti44_to_He4_Ca40_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ti44_to_He4_Ca40_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ti44_to_He4_Ca40_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_derived) = drate_dT;
        }
    }

    {
        // Ti44_to_p_Sc43_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ti44_to_p_Sc43_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ti44_to_p_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_derived) = drate_dT;
        }
    }

    {
        // He4_Ti44_to_p_V47_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ti44);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti44);
        }
#endif
        rate_He4_Ti44_to_p_V47_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ti44_to_p_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_derived) = drate_dT;
        }
    }

    {
        // Cr48_to_He4_Ti44_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cr48_to_He4_Ti44_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cr48_to_He4_Ti44_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_derived) = drate_dT;
        }
    }

    {
        // Cr48_to_p_V47_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cr48_to_p_V47_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cr48_to_p_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_derived) = drate_dT;
        }
    }

    {
        // Fe54_to_n_Fe53_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe54_to_n_Fe53_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe54_to_n_Fe53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_derived) = drate_dT;
        }
    }

    {
        // Fe53_to_n_Fe52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe53_to_n_Fe52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived) = drate_dT;
        }
    }

    {
        // Fe56_to_n_Fe55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe56_to_n_Fe55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe56_to_n_Fe55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_derived) = drate_dT;
        }
    }

    {
        // Fe55_to_n_Fe54_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe55_to_n_Fe54_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived) = drate_dT;
        }
    }

    {
        // Ni58_to_n_Ni57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni58_to_n_Ni57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni58_to_n_Ni57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_derived) = drate_dT;
        }
    }

    {
        // Ni57_to_n_Ni56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni57_to_n_Ni56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived) = drate_dT;
        }
    }


}
#endif
