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
#include <temperature_table_rates.H>
#include <tabular_weak_rates.H>
#include <derived_rates.H>

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
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_P31) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_P31) = dlog_scor_dT;
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
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ne20) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ne20) = dlog_scor_dT;
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
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 13.0_rt, 27.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Al27) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Al27) = dlog_scor_dT;
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
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 28.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Si28) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Si28) = dlog_scor_dT;
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

    reaclib_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    tabular_weak_rates::fill_rates(state.T, rhoy, Y, rate_eval);

    temp_tabular::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // fill modified rates next -- these can have ReacLib or
    // TemperatureTabular/StarLib rates as the original rate
    modified_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // fill branched rates -- these can have ReacLib or
    // TemperatureTabular/StarLib as the underlying rate or branches.
    branched_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate Derived Rates next. This should go last but before
    // approx rates, since those may have ReacLibRate,
    // TemperatureTabularRate / StarLibRate, or DerivedRate as the
    // underlying rate.
    derived_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);



    // Fill approximate rates
    approximate_rates::fill_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

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

    tabular_weak_rates::fill_rates(state.T, rhoy, Y, rate_eval);

    const auto& screened_rates = rate_eval.screened_rates;

    ydot_nuc(H1) = 0.0_rt;

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(O20) =
        (-screened_rates(k_O20_to_F20_beta_neg_weaktab)*Y(O20) + screened_rates(k_F20_to_O20_electron_capture_weaktab)*Y(F20));

    ydot_nuc(F20) =
        (-screened_rates(k_F20_to_Ne20_beta_neg_weaktab)*Y(F20) + screened_rates(k_Ne20_to_F20_electron_capture_weaktab)*Y(Ne20)) +
        (screened_rates(k_O20_to_F20_beta_neg_weaktab)*Y(O20) + -screened_rates(k_F20_to_O20_electron_capture_weaktab)*Y(F20));

    ydot_nuc(Ne20) =
        (screened_rates(k_F20_to_Ne20_beta_neg_weaktab)*Y(F20) + -screened_rates(k_Ne20_to_F20_electron_capture_weaktab)*Y(Ne20));

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Al27) = 0.0_rt;

    ydot_nuc(Si28) = 0.0_rt;

    ydot_nuc(P31) = 0.0_rt;

    ydot_nuc(S32) = 0.0_rt;

    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(H1) =
        -screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho +
        -screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho +
        0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_reaclib)*Y(Ne20)) +
        -screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho +
        -screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho +
        -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho +
        -screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho +
        0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(O16) =
        (-screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_reaclib)*Y(Ne20)) +
        -screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        -screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho;

    ydot_nuc(O20) =
        (-screened_rates(k_O20_to_F20_beta_neg_weaktab)*Y(O20) + screened_rates(k_F20_to_O20_electron_capture_weaktab)*Y(F20));

    ydot_nuc(F20) =
        (-screened_rates(k_F20_to_Ne20_beta_neg_weaktab)*Y(F20) + screened_rates(k_Ne20_to_F20_electron_capture_weaktab)*Y(Ne20)) +
        (screened_rates(k_O20_to_F20_beta_neg_weaktab)*Y(O20) + -screened_rates(k_F20_to_O20_electron_capture_weaktab)*Y(F20));

    ydot_nuc(Ne20) =
        (screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16_reaclib)*Y(Ne20)) +
        -screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho +
        (screened_rates(k_F20_to_Ne20_beta_neg_weaktab)*Y(F20) + -screened_rates(k_Ne20_to_F20_electron_capture_weaktab)*Y(Ne20));

    ydot_nuc(Mg24) =
        screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho +
        -screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho +
        (screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Al27) =
        -screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho +
        -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho +
        (-screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Si28) =
        screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho +
        screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho +
        -screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho +
        0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(P31) =
        screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho +
        -screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho +
        0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(S32) =
        screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho +
        screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho;

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
             const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates)
{

    using namespace Rates;

    amrex::Real scratch;

    scratch = -screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*state.rho - screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho - screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho;
    jac.set(H1, O16, scratch);

    jac.set(H1, O20, 0.0);

    jac.set(H1, F20, 0.0);

    jac.set(H1, Ne20, 0.0);

    scratch = screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho;
    jac.set(H1, Mg24, scratch);

    scratch = -screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho;
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho - screened_rates(k_P31_p_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(H1, P31, scratch);

    jac.set(H1, S32, 0.0);

    scratch = screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho + screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*state.rho - screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho - screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(Ne20)*state.rho - screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_Si28_He4_to_S32_reaclib)*Y(Si28)*state.rho - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = -screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho;
    jac.set(He4, O16, scratch);

    jac.set(He4, O20, 0.0);

    jac.set(He4, F20, 0.0);

    scratch = -screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16_reaclib);
    jac.set(He4, Ne20, scratch);

    scratch = -screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*state.rho - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Mg24, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(He4)*state.rho + screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = -screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*state.rho - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Si28, scratch);

    scratch = screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(He4, P31, scratch);

    jac.set(He4, S32, 0.0);

    jac.set(O16, H1, 0.0);

    scratch = -screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho - 2.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho;
    jac.set(O16, O16, scratch);

    jac.set(O16, O20, 0.0);

    jac.set(O16, F20, 0.0);

    scratch = screened_rates(k_Ne20_to_He4_O16_reaclib);
    jac.set(O16, Ne20, scratch);

    jac.set(O16, Mg24, 0.0);

    jac.set(O16, Al27, 0.0);

    jac.set(O16, Si28, 0.0);

    jac.set(O16, P31, 0.0);

    jac.set(O16, S32, 0.0);

    jac.set(O20, H1, 0.0);

    jac.set(O20, He4, 0.0);

    jac.set(O20, O16, 0.0);

    scratch = -screened_rates(k_O20_to_F20_beta_neg_weaktab);
    jac.set(O20, O20, scratch);

    scratch = screened_rates(k_F20_to_O20_electron_capture_weaktab);
    jac.set(O20, F20, scratch);

    jac.set(O20, Ne20, 0.0);

    jac.set(O20, Mg24, 0.0);

    jac.set(O20, Al27, 0.0);

    jac.set(O20, Si28, 0.0);

    jac.set(O20, P31, 0.0);

    jac.set(O20, S32, 0.0);

    jac.set(F20, H1, 0.0);

    jac.set(F20, He4, 0.0);

    jac.set(F20, O16, 0.0);

    scratch = screened_rates(k_O20_to_F20_beta_neg_weaktab);
    jac.set(F20, O20, scratch);

    scratch = -screened_rates(k_F20_to_Ne20_beta_neg_weaktab) - screened_rates(k_F20_to_O20_electron_capture_weaktab);
    jac.set(F20, F20, scratch);

    scratch = screened_rates(k_Ne20_to_F20_electron_capture_weaktab);
    jac.set(F20, Ne20, scratch);

    jac.set(F20, Mg24, 0.0);

    jac.set(F20, Al27, 0.0);

    jac.set(F20, Si28, 0.0);

    jac.set(F20, P31, 0.0);

    jac.set(F20, S32, 0.0);

    jac.set(Ne20, H1, 0.0);

    scratch = -screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(Ne20)*state.rho + screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    jac.set(Ne20, O20, 0.0);

    scratch = screened_rates(k_F20_to_Ne20_beta_neg_weaktab);
    jac.set(Ne20, F20, scratch);

    scratch = -screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*state.rho - screened_rates(k_Ne20_to_F20_electron_capture_weaktab) - screened_rates(k_Ne20_to_He4_O16_reaclib);
    jac.set(Ne20, Ne20, scratch);

    jac.set(Ne20, Mg24, 0.0);

    jac.set(Ne20, Al27, 0.0);

    jac.set(Ne20, Si28, 0.0);

    jac.set(Ne20, P31, 0.0);

    jac.set(Ne20, S32, 0.0);

    scratch = screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho + screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    jac.set(Mg24, O16, 0.0);

    jac.set(Mg24, O20, 0.0);

    jac.set(Mg24, F20, 0.0);

    scratch = screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = -screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*state.rho - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho;
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    jac.set(Mg24, Si28, 0.0);

    jac.set(Mg24, P31, 0.0);

    jac.set(Mg24, S32, 0.0);

    scratch = -screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*state.rho + screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho;
    jac.set(Al27, He4, scratch);

    jac.set(Al27, O16, 0.0);

    jac.set(Al27, O20, 0.0);

    jac.set(Al27, F20, 0.0);

    jac.set(Al27, Ne20, 0.0);

    scratch = screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(He4)*state.rho - screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Al27, Al27, scratch);

    jac.set(Al27, Si28, 0.0);

    jac.set(Al27, P31, 0.0);

    jac.set(Al27, S32, 0.0);

    scratch = screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*state.rho + screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_Si28_He4_to_S32_reaclib)*Y(Si28)*state.rho - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    jac.set(Si28, O20, 0.0);

    jac.set(Si28, F20, 0.0);

    jac.set(Si28, Ne20, 0.0);

    scratch = screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = -screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*state.rho - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho;
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    jac.set(Si28, S32, 0.0);

    scratch = -screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho - screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*state.rho + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    jac.set(P31, O20, 0.0);

    jac.set(P31, F20, 0.0);

    jac.set(P31, Ne20, 0.0);

    jac.set(P31, Mg24, 0.0);

    scratch = screened_rates(k_Al27_He4_to_P31_reaclib)*Y(He4)*state.rho;
    jac.set(P31, Al27, scratch);

    scratch = screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    scratch = -screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho - screened_rates(k_P31_p_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(P31, P31, scratch);

    jac.set(P31, S32, 0.0);

    scratch = screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = screened_rates(k_Si28_He4_to_S32_reaclib)*Y(Si28)*state.rho;
    jac.set(S32, He4, scratch);

    jac.set(S32, O16, 0.0);

    jac.set(S32, O20, 0.0);

    jac.set(S32, F20, 0.0);

    jac.set(S32, Ne20, 0.0);

    jac.set(S32, Mg24, 0.0);

    jac.set(S32, Al27, 0.0);

    scratch = screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = screened_rates(k_P31_p_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);

    jac.set(S32, S32, 0.0);


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
