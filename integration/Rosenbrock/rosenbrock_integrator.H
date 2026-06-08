#ifndef ROSENBROCK_INTEGRATOR_H
#define ROSENBROCK_INTEGRATOR_H

#include <AMReX_Algorithm.H>

#include <cmath>
#include <limits>

#include <network.H>
#include <burn_type.H>
#include <linpack.H>

#ifdef STRANG
#include <integrator_rhs_strang.H>
#endif
#ifdef SDC
#include <integrator_rhs_sdc.H>
#endif
#ifdef NSE_TABLE
#include <nse_table_check.H>
#endif
#ifdef NSE_NET
#include <nse_check.H>
#endif

#include <integrator_data.H>
#include <rosenbrock_type.H>
#include <rosenbrock_tableau.H>

namespace rosenbrock {

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real& stage (rosenbrock_t<int_neqs>& rstate, const int i, const int n)
{
    if (i == 1) {
        return rstate.ak1(n);
    }
    if (i == 2) {
        return rstate.ak2(n);
    }
    if (i == 3) {
        return rstate.ak3(n);
    }
    if (i == 4) {
        return rstate.ak4(n);
    }
    if (i == 5) {
        return rstate.ak5(n);
    }
    if (i == 6) {
        return rstate.ak6(n);
    }
    if (i == 7) {
        return rstate.ak7(n);
    }
    if (i == 8) {
        return rstate.ak8(n);
    }
    return rstate.work(n);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real stage (const rosenbrock_t<int_neqs>& rstate, const int i, const int n)
{
    if (i == 1) {
        return rstate.ak1(n);
    }
    if (i == 2) {
        return rstate.ak2(n);
    }
    if (i == 3) {
        return rstate.ak3(n);
    }
    if (i == 4) {
        return rstate.ak4(n);
    }
    if (i == 5) {
        return rstate.ak5(n);
    }
    if (i == 6) {
        return rstate.ak6(n);
    }
    if (i == 7) {
        return rstate.ak7(n);
    }
    if (i == 8) {
        return rstate.ak8(n);
    }
    return rstate.work(n);
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RArray1D& stage_vector (rosenbrock_t<int_neqs>& rstate, const int i)
{
    if (i == 1) {
        return rstate.ak1;
    }
    if (i == 2) {
        return rstate.ak2;
    }
    if (i == 3) {
        return rstate.ak3;
    }
    if (i == 4) {
        return rstate.ak4;
    }
    if (i == 5) {
        return rstate.ak5;
    }
    if (i == 6) {
        return rstate.ak6;
    }
    if (i == 7) {
        return rstate.ak7;
    }
    if (i == 8) {
        return rstate.ak8;
    }
    return rstate.work;
}

template <typename RosenbrockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rtol_for (const RosenbrockT& rstate, const int n)
{
    return n == net_ienuc ? rstate.rtol_enuc : rstate.rtol_spec;
}

template <typename RosenbrockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real atol_for (const RosenbrockT& rstate, const int n)
{
    return n == net_ienuc ? rstate.atol_enuc : rstate.atol_spec;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real error_norm (const rosenbrock_t<int_neqs>& rstate)
{
    amrex::Real err = 0.0_rt;
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real sk = atol_for(rstate, n) +
            rtol_for(rstate, n) * amrex::max(std::abs(rstate.y(n)), std::abs(rstate.ynew(n)));
        const amrex::Real term = rstate.work(n) / sk;
        err += term * term;
    }
    return std::sqrt(err / static_cast<amrex::Real>(int_neqs));
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool species_change_valid (const rosenbrock_t<int_neqs>& rstate)
{
#ifdef STRANG
    constexpr amrex::Real increase_change_factor = 4.0_rt;
    constexpr amrex::Real decrease_change_factor = 0.25_rt;

    for (int i = 1; i <= NumSpec; ++i) {
        if (std::abs(rstate.y(i)) > integrator_rp::X_reject_buffer * rstate.atol_spec &&
            std::abs(rstate.ynew(i)) > integrator_rp::X_reject_buffer * rstate.atol_spec &&
            (std::abs(rstate.ynew(i)) > increase_change_factor * std::abs(rstate.y(i)) ||
             std::abs(rstate.ynew(i)) < decrease_change_factor * std::abs(rstate.y(i)))) {
            return false;
        }
    }
#endif

    return true;
}

template <typename MatrixType, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool matrix_is_finite (const MatrixType& matrix)
{
    for (int j = 1; j <= int_neqs; ++j) {
        for (int i = 1; i <= int_neqs; ++i) {
            const amrex::Real value = matrix(i, j);
            if (std::isnan(value) || std::isinf(value)) {
                return false;
            }
        }
    }

    return true;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_integrator_state (const amrex::Array1D<amrex::Real, 1, int_neqs>& y)
{
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real yn = y(n);
        if (std::isnan(yn) || std::isinf(yn)) {
            return false;
        }
    }

#ifdef STRANG
    for (int n = 1; n <= NumSpec; ++n) {
        if (y(n) < -integrator_rp::atol_spec) {
            return false;
        }

        if (! integrator_rp::use_number_densities && y(n) > 1.0_rt) {
            return false;
        }
    }

    if (y(net_ienuc) <= 0.0_rt) {
        return false;
    }
#endif

    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_integrator_component (const int n, const amrex::Real yn)
{
    if (std::isnan(yn) || std::isinf(yn)) {
        return false;
    }

#ifdef STRANG
    if (n <= NumSpec) {
        if (yn < -integrator_rp::atol_spec) {
            return false;
        }

        if (! integrator_rp::use_number_densities && yn > 1.0_rt) {
            return false;
        }
    }

    if (n == net_ienuc && yn <= 0.0_rt) {
        return false;
    }
#endif

    return true;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_rodas5p_stage5 (const rodas5p_t<int_neqs>& rstate)
{
    using C = rodas5p_tableau;

    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real yn = rstate.y(n) +
            C::a(5, 1) * rstate.ak1(n) +
            C::a(5, 2) * rstate.ak2(n) +
            C::a(5, 3) * rstate.ak3(n) +
            C::a(5, 4) * rstate.ak4(n);
        if (! valid_integrator_component(n, yn)) {
            return false;
        }
    }

    return true;
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_rodas5p_update_stage (const rodas5p_t<int_neqs>& rstate,
                                 const amrex::Array1D<amrex::Real, 1, int_neqs>& update)
{
    for (int n = 1; n <= int_neqs; ++n) {
        if (! valid_integrator_component(n, rstate.y(n) + update(n))) {
            return false;
        }
    }

    return true;
}

#ifdef STRANG
template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void finish_rodas5p_stage_state (const amrex::Real time, BurnT& state)
{
    if (integrator_rp::do_species_clip) {
        for (int n = 1; n <= NumSpec; ++n) {
            state.xn[n-1] = amrex::Clamp(state.xn[n-1], integrator_rp::SMALL_X_SAFE, 1.0_rt);
        }
    }

    if (integrator_rp::use_number_densities) {
        for (int n = 1; n <= NumSpec; ++n) {
            state.xn[n-1] = amrex::max(state.xn[n-1], integrator_rp::SMALL_X_SAFE);
        }
    }

    if (integrator_rp::renormalize_abundances) {
        amrex::Real sum = 0.0_rt;

        for (int n = 1; n <= NumSpec; ++n) {
            sum += state.xn[n-1];
        }

        for (int n = 1; n <= NumSpec; ++n) {
            state.xn[n-1] /= sum;
        }
    }

#ifdef AUX_THERMO
    state.abar = 0.0_rt;
    state.y_e = 0.0_rt;
    for (int n = 1; n <= NumSpec; ++n) {
      state.abar += state.xn[n-1] * aion_inv[n-1];
      state.y_e += state.xn[n-1] * zion[n-1] * aion_inv[n-1];
    }
    state.abar = 1.0_rt / state.abar;
    state.zbar = state.abar * state.y_e;

    state.aux[AuxZero::iabar] = state.abar;
    state.aux[AuxZero::iye] = state.y_e;
    state.aux[AuxZero::ibea] = 0.0_rt;
#endif

    if (integrator_rp::scale_system) {
        state.e *= state.e_scale;
    }

    if (integrator_rp::call_eos_in_rhs) {
        eos(eos_input_re, state);
    }

    if (state.T_fixed > 0.0_rt) {
        state.T = state.T_fixed;
    }

    state.time = time;
}
#endif

#ifdef SDC
template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void finish_rodas5p_stage_state (const amrex::Real time, BurnT& state)
{
    update_density_in_time(time, state);

    if (integrator_rp::do_species_clip) {
        for (int n = 1; n <= NumSpec; ++n) {
            state.y[SFS+n-1] = amrex::Clamp(state.y[SFS+n-1],
                                            state.rho * integrator_rp::SMALL_X_SAFE,
                                            state.rho);
        }
    }

    if (integrator_rp::renormalize_abundances) {
        amrex::Real nspec_sum = 0.0_rt;
        for (int n = 1; n <= NumSpec; n++) {
            nspec_sum += state.y[SFS+n-1];
        }
        nspec_sum /= state.y[SRHO];

        for (int n = 1; n <= NumSpec; n++) {
            state.y[SFS+n-1] /= nspec_sum;
        }
    }

    if (integrator_rp::scale_system) {
        state.y[SEINT] *= state.e_scale;
    }

    update_density_in_time(time, state);

    amrex::Real rhoInv = 1.0_rt / state.rho;

    for (int n = 0; n < NumSpec; n++) {
        state.xn[n] = state.y[SFS+n] * rhoInv;
    }

#ifdef AUX_THERMO
    set_aux_comp_from_X(state);

    for (int n = 0; n < NumAux; n++) {
        state.y[SFX+n] = state.rho * state.aux[n];
    }
#endif

    state.e = state.y[SEINT] * rhoInv;

    if (integrator_rp::call_eos_in_rhs) {
        eos(eos_input_re, state);
    }

    if (state.T_fixed > 0.0_rt) {
        state.T = state.T_fixed;
    }

    state.time = time;
}
#endif

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void prepare_rodas5p_stage5_state (const amrex::Real time, BurnT& state,
                                   const rodas5p_t<int_neqs>& rstate)
{
    using C = rodas5p_tableau;

#ifdef STRANG
    for (int n = 1; n <= NumSpec; ++n) {
        const amrex::Real yn = rstate.y(n) +
            C::a(5, 1) * rstate.ak1(n) +
            C::a(5, 2) * rstate.ak2(n) +
            C::a(5, 3) * rstate.ak3(n) +
            C::a(5, 4) * rstate.ak4(n);
        state.y[n-1] = yn;
        state.xn[n-1] = yn;
    }

    state.y[SEINT] = rstate.y(net_ienuc) +
        C::a(5, 1) * rstate.ak1(net_ienuc) +
        C::a(5, 2) * rstate.ak2(net_ienuc) +
        C::a(5, 3) * rstate.ak3(net_ienuc) +
        C::a(5, 4) * rstate.ak4(net_ienuc);
    state.e = state.y[SEINT];
#endif

#ifdef SDC
    for (int n = 0; n < SVAR_EVOLVE; n++) {
        const int i = n + 1;
        state.y[n] = rstate.y(i) +
            C::a(5, 1) * rstate.ak1(i) +
            C::a(5, 2) * rstate.ak2(i) +
            C::a(5, 3) * rstate.ak3(i) +
            C::a(5, 4) * rstate.ak4(i);
    }
#endif

    finish_rodas5p_stage_state(time, state);
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void prepare_rodas5p_update_stage_state (const amrex::Real time, BurnT& state,
                                         const rodas5p_t<int_neqs>& rstate,
                                         const amrex::Array1D<amrex::Real, 1, int_neqs>& update)
{
#ifdef STRANG
    for (int n = 1; n <= NumSpec; ++n) {
        const amrex::Real yn = rstate.y(n) + update(n);
        state.y[n-1] = yn;
        state.xn[n-1] = yn;
    }

    state.y[SEINT] = rstate.y(net_ienuc) + update(net_ienuc);
    state.e = state.y[SEINT];
#endif

#ifdef SDC
    for (int n = 0; n < SVAR_EVOLVE; n++) {
        const int i = n + 1;
        state.y[n] = rstate.y(i) + update(i);
    }
#endif

    finish_rodas5p_stage_state(time, state);
}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rodas5p_stage_component (const BurnT& state, const int n)
{
#ifdef STRANG
    if (n <= NumSpec) {
        return state.y[n-1];
    }

    return state.y[SEINT];
#endif

#ifdef SDC
    const int i = n - 1;
    if (i == SEINT && integrator_rp::scale_system) {
        return state.y[i] / state.e_scale;
    }

    return state.y[i];
#endif
}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rodas5p_rhs_from_prepared_state (const amrex::Real time, BurnT& state,
                                      RArray1D& ydot,
                                      [[maybe_unused]] const bool in_jacobian=false)
{
    amrex::ignore_unused(time);

    if (state.T <= EOSData::mintemp || state.T >= integrator_rp::MAX_TEMP) {
        for (int n = 1; n <= INT_NEQS; ++n) {
            ydot(n) = 0.0_rt;
        }

        return;
    }

#ifdef STRANG
#ifdef NEW_NETWORK_IMPLEMENTATION
    RHS::rhs(state, ydot);
#else
    actual_rhs(state, ydot);
#endif

#ifdef NONAKA_PLOT
    if (! in_jacobian) {
        nonaka_rhs(time, state, ydot);
    }
#endif

    if (! integrator_rp::use_number_densities) {
        for (int n = 1; n <= NumSpec; ++n) {
            ydot(n) *= aion[n-1];
        }
    }

    if (integrator_rp::scale_system) {
        ydot(net_ienuc) /= state.e_scale;
    }

    if (! integrator_rp::integrate_energy) {
        ydot(net_ienuc) = 0.0_rt;
    }

    if (integrator_rp::react_boost > 0.0_rt) {
        for (int n = 1; n <= INT_NEQS; ++n) {
            ydot(n) *= integrator_rp::react_boost;
        }
    }
#endif

#ifdef SDC
    actual_rhs(state, ydot);

#ifdef NONAKA_PLOT
    if (! in_jacobian) {
        nonaka_rhs(time, state, ydot);
    }
#endif

    if (integrator_rp::react_boost > 0.0_rt) {
        for (int n = 1; n <= neqs; ++n) {
            ydot(n) *= integrator_rp::react_boost;
        }
    }

    rhs_to_int(time, state, ydot);
#endif
}

template <int int_neqs>
struct rhs_array_ref {
    amrex::Array1D<amrex::Real, 1, int_neqs>* data;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real& operator() (const int n) const {
        return (*data)(n);
    }
};

template <int int_neqs>
struct rhs_state_view {
    amrex::Real t;
    rhs_array_ref<int_neqs> y;
};

template <typename BurnT, typename RosenbrockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs_at_array (const amrex::Real time, BurnT& state, RosenbrockT& rstate,
                   amrex::Array1D<amrex::Real, 1, RosenbrockT::neqs>& y,
                   RArray1D& ydot)
{
    amrex::ignore_unused(rstate);
    constexpr int int_neqs = RosenbrockT::neqs;
    rhs_state_view<int_neqs> rhs_state{time, rhs_array_ref<int_neqs>{&y}};

    rhs(time, state, rhs_state, ydot);
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_jacobian (BurnT& state, rosenbrock_t<int_neqs>& rstate, const amrex::Real time)
{
    constexpr bool in_jacobian = true;

    rhs(time, state, rstate, rstate.ak1, in_jacobian);
    rstate.n_rhs += 1;

    if (rstate.jacobian_type == 1) {
        jac(time, state, rstate, rstate.jac);
        rstate.n_jac += 1;

        if (matrix_is_finite<decltype(rstate.jac), int_neqs>(rstate.jac)) {
            return;
        }
    }

    constexpr amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Array1D<amrex::Real, 1, int_neqs> ewt;
    amrex::Array1D<amrex::Real, 1, int_neqs> ybase;

    // rhs() may clean or renormalize all components of rstate.y.
    for (int i = 1; i <= int_neqs; ++i) {
        ybase(i) = rstate.y(i);
    }

    amrex::Real fac = 0.0_rt;
    for (int i = 1; i <= int_neqs; ++i) {
        ewt(i) = 1.0_rt / (rtol_for(rstate, i) * std::abs(rstate.y(i)) + atol_for(rstate, i));
        fac += (rstate.ak1(i) * ewt(i)) * (rstate.ak1(i) * ewt(i));
    }
    fac = std::sqrt(fac / static_cast<amrex::Real>(int_neqs));

    amrex::Real R0 = 1000.0_rt * std::abs(rstate.dt) *
        UROUND * static_cast<amrex::Real>(int_neqs) * fac;
    if (R0 == 0.0_rt) {
        R0 = 1.0_rt;
    }

    for (int j = 1; j <= int_neqs; ++j) {
        const amrex::Real yj = ybase(j);
        const amrex::Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / ewt(j));
        for (int n = 1; n <= int_neqs; ++n) {
            rstate.y(n) = ybase(n);
        }
        rstate.y(j) += R;

        rhs(time, state, rstate, rstate.work, in_jacobian);

        fac = 1.0_rt / R;
        for (int i = 1; i <= int_neqs; ++i) {
            rstate.jac.set(i, j, (rstate.work(i) - rstate.ak1(i)) * fac);
        }

        for (int n = 1; n <= int_neqs; ++n) {
            rstate.y(n) = ybase(n);
        }
    }

    rstate.n_rhs += int_neqs;
    rstate.n_jac += 1;
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_jacobian_rodas5p (BurnT& state, rodas5p_t<int_neqs>& rstate,
                                const amrex::Real time)
{
    constexpr bool in_jacobian = true;

    rhs(time, state, rstate, rstate.ak1, in_jacobian);
    rstate.n_rhs += 1;

    if (rstate.jacobian_type == 1) {
        jac(time, state, rstate, rstate.jac);
        rstate.n_jac += 1;

        if (matrix_is_finite<decltype(rstate.jac), int_neqs>(rstate.jac)) {
            return;
        }
    }

    constexpr amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();

    // rhs() may clean or renormalize all components of rstate.y.
    for (int i = 1; i <= int_neqs; ++i) {
        rstate.ak2(i) = rstate.y(i);
    }

    amrex::Real fac = 0.0_rt;
    for (int i = 1; i <= int_neqs; ++i) {
        const amrex::Real ewt =
            1.0_rt / (rtol_for(rstate, i) * std::abs(rstate.y(i)) + atol_for(rstate, i));
        fac += (rstate.ak1(i) * ewt) * (rstate.ak1(i) * ewt);
    }
    fac = std::sqrt(fac / static_cast<amrex::Real>(int_neqs));

    amrex::Real R0 = 1000.0_rt * std::abs(rstate.dt) *
        UROUND * static_cast<amrex::Real>(int_neqs) * fac;
    if (R0 == 0.0_rt) {
        R0 = 1.0_rt;
    }

    for (int j = 1; j <= int_neqs; ++j) {
        const amrex::Real yj = rstate.ak2(j);
        const amrex::Real ewt_j =
            1.0_rt / (rtol_for(rstate, j) * std::abs(yj) + atol_for(rstate, j));
        const amrex::Real R = amrex::max(std::sqrt(UROUND) * std::abs(yj), R0 / ewt_j);
        for (int n = 1; n <= int_neqs; ++n) {
            rstate.y(n) = rstate.ak2(n);
        }
        rstate.y(j) += R;

        rhs(time, state, rstate, rstate.ak3, in_jacobian);

        fac = 1.0_rt / R;
        for (int i = 1; i <= int_neqs; ++i) {
            rstate.jac.set(i, j, (rstate.ak3(i) - rstate.ak1(i)) * fac);
        }
    }

    for (int n = 1; n <= int_neqs; ++n) {
        rstate.y(n) = rstate.ak2(n);
    }

    rstate.n_rhs += int_neqs;
    rstate.n_jac += 1;
}

template <typename BurnT, int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs_at (const amrex::Real time, BurnT& state, rosenbrock_t<int_neqs>& rstate,
             const amrex::Array1D<amrex::Real, 1, int_neqs>& y,
             amrex::Array1D<amrex::Real, 1, int_neqs>& ydot)
{
    amrex::Array1D<amrex::Real, 1, int_neqs> ysave;

    for (int n = 1; n <= int_neqs; ++n) {
        ysave(n) = rstate.y(n);
        rstate.y(n) = y(n);
    }

    rhs(time, state, rstate, ydot);

    for (int n = 1; n <= int_neqs; ++n) {
        rstate.y(n) = ysave(n);
    }
}

template <typename RosenbrockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int decompose (RosenbrockT& rstate, const amrex::Real fac)
{
    constexpr int int_neqs = RosenbrockT::neqs;

    for (int j = 1; j <= int_neqs; ++j) {
        for (int i = 1; i <= int_neqs; ++i) {
            rstate.jac(i, j) = -rstate.jac(i, j);
        }
        rstate.jac(j, j) += fac;
    }

    int ierr_linpack = 0;
    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgefa<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, ierr_linpack);
    } else {
        constexpr bool allow_pivot{false};
        dgefa<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, ierr_linpack);
    }

    return ierr_linpack;
}

template <typename RosenbrockT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve (RosenbrockT& rstate, RArray1D& b)
{
    constexpr int int_neqs = RosenbrockT::neqs;

    if (integrator_rp::linalg_do_pivoting == 1) {
        constexpr bool allow_pivot{true};
        dgesl<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, b);
    } else {
        constexpr bool allow_pivot{false};
        dgesl<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, b);
    }
}

} // namespace rosenbrock

template <typename BurnT, typename Tableau>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rosenbrock_integrator (BurnT& state, rosenbrock_t<integrator_neqs<BurnT>()>& rstate)
{
    using C = Tableau;
    constexpr int int_neqs = integrator_neqs<BurnT>();
    static_assert(C::stages <= 8, "Rosenbrock integrator stores up to 8 stages");

    if (rstate.tout == rstate.t) {
        return IERR_SUCCESS;
    }

    for (int n = 1; n <= int_neqs; ++n) {
        if (rosenbrock::atol_for(rstate, n) <= 0.0_rt ||
            rosenbrock::rtol_for(rstate, n) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
            return IERR_TOO_MUCH_ACCURACY_REQUESTED;
        }
    }

    rstate.n_rhs = 0;
    rstate.n_jac = 0;
    rstate.n_step = 0;

    rstate.dt = rstate.tout - rstate.t;

    const amrex::Real posneg = rstate.tout >= rstate.t ? 1.0_rt : -1.0_rt;
    const amrex::Real hmax = amrex::min(integrator_rp::ode_max_dt,
                                        std::abs(rstate.tout - rstate.t));
    amrex::Real h = amrex::min(std::abs(rstate.dt), hmax) * posneg;
    if (std::abs(h) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
        h = 1.e-6_rt * posneg;
    }

    bool reject = false;
    bool last = false;
    int nsing = 0;
    int n_reject = 0;
    amrex::Real facold = 1.0_rt;
    amrex::Real errold = 1.0_rt;
    amrex::Real hopt = h;
    amrex::Real x = rstate.t;
    amrex::Array1D<amrex::Real, 1, int_neqs> rhs_tmp;

    // main evolution loop

    // within this loop, x will be the current time and h will be the
    // step size we are attempting.  The solution at the current time
    // is rstate.y().

    while (true) {

        if (rstate.n_step > integrator_rp::ode_max_steps) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_TOO_MANY_STEPS;
        }

        if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_DT_UNDERFLOW;
        }

        if (last) {
            rstate.t = x;
            rstate.dt = hopt;
            return IERR_SUCCESS;
        }

        hopt = h;
        if ((x + h * (1.0_rt + timestep_safety_factor) - rstate.tout) * posneg >= 0.0_rt) {
            h = rstate.tout - x;
            last = true;
        }

#ifdef NSE
        // check if, during the course of integration, we hit NSE, and
        // if so, bail out.

        // we only do this after MIN_NSE_BAILOUT_STEPS to prevent us
        // from hitting this right at the start when the integrator
        // might do so wild exploration.  Also ensure we are not
        // working > tmax, so we don't need to worry about
        // extrapolating back in time.

        if (rstate.n_step > MIN_NSE_BAILOUT_STEPS && x <= rstate.tout) {
            // first we need to make the burn_t in sync

            int_to_burn(x, rstate, state);

            if (in_nse(state)) {
                rstate.t = x;
                return IERR_ENTERED_NSE;
            }
        }
#endif


        // take a step -- this loop will keep trying until a step is
        // successful or a catastrophic error is encountered.

        while (true) {
            rstate.dt = h;

            if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
                rstate.t = x;
                rstate.dt = h;
                return IERR_DT_UNDERFLOW;
            }

            const amrex::Real fac = 1.0_rt / (h * C::gamma);
            rosenbrock::evaluate_jacobian(state, rstate, x);
            int ierr_linpack = rosenbrock::decompose(rstate, fac);

            if (ierr_linpack != 0) {
                nsing += 1;
                if (nsing >= 5) {
                    rstate.t = x;
                    rstate.dt = h;
                    return IERR_LU_DECOMPOSITION_ERROR;
                }
                h *= 0.5_rt;
                reject = true;
                last = false;
                continue;
            }

            rosenbrock::solve(rstate, rstate.ak1);

            bool valid_stage_state = true;
            for (int istage = 2; istage <= C::stages; ++istage) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ynew(n) = rstate.y(n);
                    rosenbrock::stage(rstate, istage, n) = 0.0_rt;

                    for (int jstage = 1; jstage < istage; ++jstage) {
                        const amrex::Real akj = rosenbrock::stage(rstate, jstage, n);
                        rstate.ynew(n) += C::a(istage, jstage) * akj;
                        rosenbrock::stage(rstate, istage, n) +=
                            (C::c(istage, jstage) / h) * akj;
                    }
                }

                if (! rosenbrock::valid_integrator_state(rstate.ynew)) {
                    h *= 0.25_rt;
                    reject = true;
                    n_reject += 1;
                    last = false;
                    valid_stage_state = false;
                    break;
                }

                rosenbrock::rhs_at(x + C::ct(istage) * h, state, rstate, rstate.ynew, rhs_tmp);
                rstate.n_rhs += 1;

                for (int n = 1; n <= int_neqs; ++n) {
                    rosenbrock::stage(rstate, istage, n) += rhs_tmp(n);
                }
                rosenbrock::solve(rstate, rosenbrock::stage_vector(rstate, istage));
            }

            if (! valid_stage_state) {
                continue;
            }

            for (int n = 1; n <= int_neqs; ++n) {
                amrex::Real y_update = 0.0_rt;
                amrex::Real err = 0.0_rt;

                for (int istage = 1; istage <= C::stages; ++istage) {
                    const amrex::Real aki = rosenbrock::stage(rstate, istage, n);
                    y_update += C::b(istage) * aki;
                    err += C::e(istage) * aki;
                }

                rstate.ynew(n) = rstate.y(n) + y_update;
                rstate.work(n) = err;
            }

            rstate.n_step += 1;

            const bool valid_state = rosenbrock::valid_integrator_state(rstate.ynew);
            if (! valid_state) {
                reject = true;
                n_reject += 1;
                last = false;
                h *= 0.25_rt;
                continue;
            }

            const bool valid_update = rosenbrock::species_change_valid(rstate);
            constexpr amrex::Real err_min = 1.e-10_rt;
            amrex::Real err = rosenbrock::error_norm(rstate);
            amrex::Real controller_fac = amrex::Clamp(
                std::pow(1.0_rt / amrex::max(err, err_min),
                         1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                std::pow(1.0_rt / amrex::max(errold, err_min),
                         1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                std::pow(facold, -1.0_rt / integrator_rp::h211b_b),
                integrator_rp::h211b_fac_min, integrator_rp::h211b_fac_max);
            facold = controller_fac;
            errold = err;

            amrex::Real hnew = h * controller_fac;

            // update the starting time and state for the next step
            // and predict the next timestep

            if (err <= 1.0_rt && valid_update) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.y(n) = rstate.ynew(n);
                }
                x += h;

                if (std::abs(hnew) > hmax) {
                    hnew = posneg * hmax;
                }
                if (reject) {
                    hnew = posneg * amrex::min(std::abs(hnew), std::abs(h));
                }
                reject = false;
                n_reject = 0;
                h = hnew;
                break;
            }

            // if we made it here, then we've failed
            // cut the timestep in prep for trying again

            reject = true;
            n_reject += 1;
            last = false;
            if (n_reject >= 2) {
                hnew *= integrator_rp::h211b_reduction_fac;
            }
            if (valid_update) {
                h = posneg * amrex::min(std::abs(hnew),
                                         integrator_rp::h211b_reduction_fac * std::abs(h));
            } else {
                h *= 0.25_rt;
            }
        }
    }
}

namespace rosenbrock {

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real error_norm (const rodas5p_t<int_neqs>& rstate,
                        const amrex::Array1D<amrex::Real, 1, int_neqs>& ynew,
                        const amrex::Array1D<amrex::Real, 1, int_neqs>& error)
{
    amrex::Real err = 0.0_rt;
    for (int n = 1; n <= int_neqs; ++n) {
        const amrex::Real sk = atol_for(rstate, n) +
            rtol_for(rstate, n) * amrex::max(std::abs(rstate.y(n)), std::abs(ynew(n)));
        const amrex::Real term = error(n) / sk;
        err += term * term;
    }
    return std::sqrt(err / static_cast<amrex::Real>(int_neqs));
}

template <int int_neqs>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool species_change_valid (const rodas5p_t<int_neqs>& rstate,
                           const amrex::Array1D<amrex::Real, 1, int_neqs>& ynew)
{
#ifdef STRANG
    constexpr amrex::Real increase_change_factor = 4.0_rt;
    constexpr amrex::Real decrease_change_factor = 0.25_rt;

    for (int i = 1; i <= NumSpec; ++i) {
        if (std::abs(rstate.y(i)) > integrator_rp::X_reject_buffer * rstate.atol_spec &&
            std::abs(ynew(i)) > integrator_rp::X_reject_buffer * rstate.atol_spec &&
            (std::abs(ynew(i)) > increase_change_factor * std::abs(rstate.y(i)) ||
             std::abs(ynew(i)) < decrease_change_factor * std::abs(rstate.y(i)))) {
            return false;
        }
    }
#endif

    return true;
}

} // namespace rosenbrock

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rodas5p_integrator (BurnT& state, rodas5p_t<integrator_neqs<BurnT>()>& rstate)
{
    using C = rosenbrock::rodas5p_tableau;
    constexpr int int_neqs = integrator_neqs<BurnT>();

    if (rstate.tout == rstate.t) {
        return IERR_SUCCESS;
    }

    for (int n = 1; n <= int_neqs; ++n) {
        if (rosenbrock::atol_for(rstate, n) <= 0.0_rt ||
            rosenbrock::rtol_for(rstate, n) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
            return IERR_TOO_MUCH_ACCURACY_REQUESTED;
        }
    }

    rstate.n_rhs = 0;
    rstate.n_jac = 0;
    rstate.n_step = 0;

    rstate.dt = rstate.tout - rstate.t;

    const amrex::Real posneg = rstate.tout >= rstate.t ? 1.0_rt : -1.0_rt;
    const amrex::Real hmax = amrex::min(integrator_rp::ode_max_dt,
                                        std::abs(rstate.tout - rstate.t));
    amrex::Real h = amrex::min(std::abs(rstate.dt), hmax) * posneg;
    if (std::abs(h) <= 10.0_rt * std::numeric_limits<amrex::Real>::epsilon()) {
        h = 1.e-6_rt * posneg;
    }

    bool reject = false;
    bool last = false;
    int nsing = 0;
    int n_reject = 0;
    amrex::Real facold = 1.0_rt;
    amrex::Real errold = 1.0_rt;
    amrex::Real hopt = h;
    amrex::Real x = rstate.t;

    while (true) {

        if (rstate.n_step > integrator_rp::ode_max_steps) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_TOO_MANY_STEPS;
        }

        if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
            rstate.t = x;
            rstate.dt = h;
            return IERR_DT_UNDERFLOW;
        }

        if (last) {
            rstate.t = x;
            rstate.dt = hopt;
            return IERR_SUCCESS;
        }

        hopt = h;
        if ((x + h * (1.0_rt + timestep_safety_factor) - rstate.tout) * posneg >= 0.0_rt) {
            h = rstate.tout - x;
            last = true;
        }

#ifdef NSE
        if (rstate.n_step > MIN_NSE_BAILOUT_STEPS && x <= rstate.tout) {
            int_to_burn(x, rstate, state);

            if (in_nse(state)) {
                rstate.t = x;
                return IERR_ENTERED_NSE;
            }
        }
#endif

        while (true) {
            rstate.dt = h;

            if (0.1_rt * std::abs(h) <= std::abs(x) * std::numeric_limits<amrex::Real>::epsilon()) {
                rstate.t = x;
                rstate.dt = h;
                return IERR_DT_UNDERFLOW;
            }

            const amrex::Real fac = 1.0_rt / (h * C::gamma);
            rosenbrock::evaluate_jacobian_rodas5p(state, rstate, x);
            int ierr_linpack = rosenbrock::decompose(rstate, fac);

            if (ierr_linpack != 0) {
                nsing += 1;
                if (nsing >= 5) {
                    rstate.t = x;
                    rstate.dt = h;
                    return IERR_LU_DECOMPOSITION_ERROR;
                }
                h *= 0.5_rt;
                reject = true;
                last = false;
                continue;
            }

            rosenbrock::solve(rstate, rstate.ak1);

            bool valid_stage_state = true;

            for (int n = 1; n <= int_neqs; ++n) {
                rstate.ak2(n) = rstate.y(n) + C::a(2, 1) * rstate.ak1(n);
            }
            valid_stage_state = rosenbrock::valid_integrator_state(rstate.ak2);
            if (valid_stage_state) {
                rosenbrock::rhs_at_array(x + C::ct(2) * h, state, rstate, rstate.ak2, rstate.ak2);
                rstate.n_rhs += 1;
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak2(n) += (C::c(2, 1) / h) * rstate.ak1(n);
                }
                rosenbrock::solve(rstate, rstate.ak2);
            }

            if (valid_stage_state) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak3(n) = rstate.y(n) +
                        C::a(3, 1) * rstate.ak1(n) +
                        C::a(3, 2) * rstate.ak2(n);
                }
                valid_stage_state = rosenbrock::valid_integrator_state(rstate.ak3);
            }
            if (valid_stage_state) {
                rosenbrock::rhs_at_array(x + C::ct(3) * h, state, rstate, rstate.ak3, rstate.ak3);
                rstate.n_rhs += 1;
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak3(n) +=
                        (C::c(3, 1) * rstate.ak1(n) +
                         C::c(3, 2) * rstate.ak2(n)) / h;
                }
                rosenbrock::solve(rstate, rstate.ak3);
            }

            if (valid_stage_state) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak4(n) = rstate.y(n) +
                        C::a(4, 1) * rstate.ak1(n) +
                        C::a(4, 2) * rstate.ak2(n) +
                        C::a(4, 3) * rstate.ak3(n);
                }
                valid_stage_state = rosenbrock::valid_integrator_state(rstate.ak4);
            }
            if (valid_stage_state) {
                rosenbrock::rhs_at_array(x + C::ct(4) * h, state, rstate, rstate.ak4, rstate.ak4);
                rstate.n_rhs += 1;
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak4(n) +=
                        (C::c(4, 1) * rstate.ak1(n) +
                         C::c(4, 2) * rstate.ak2(n) +
                         C::c(4, 3) * rstate.ak3(n)) / h;
                }
                rosenbrock::solve(rstate, rstate.ak4);
            }

            if (valid_stage_state) {
                valid_stage_state = rosenbrock::valid_rodas5p_stage5(rstate);
            }
            if (valid_stage_state) {
                rosenbrock::prepare_rodas5p_stage5_state(x + C::ct(5) * h, state, rstate);
                rosenbrock::rodas5p_rhs_from_prepared_state(x + C::ct(5) * h, state, rstate.ak4);
                rstate.n_rhs += 1;
                for (int n = 1; n <= int_neqs; ++n) {
                    const amrex::Real k4 =
                        (rosenbrock::rodas5p_stage_component(state, n) -
                         rstate.y(n) -
                         C::a(5, 1) * rstate.ak1(n) -
                         C::a(5, 2) * rstate.ak2(n) -
                         C::a(5, 3) * rstate.ak3(n)) / C::a(5, 4);

                    rstate.ak4(n) +=
                        (C::c(5, 1) * rstate.ak1(n) +
                         C::c(5, 2) * rstate.ak2(n) +
                         C::c(5, 3) * rstate.ak3(n) +
                         C::c(5, 4) * k4) / h;
                }
                rosenbrock::solve(rstate, rstate.ak4);
            }

            if (! valid_stage_state) {
                h *= 0.25_rt;
                reject = true;
                n_reject += 1;
                last = false;
                continue;
            }

            for (int n = 1; n <= int_neqs; ++n) {
                const amrex::Real k1 = rstate.ak1(n);
                const amrex::Real k2 = rstate.ak2(n);
                const amrex::Real k3 = rstate.ak3(n);
                const amrex::Real k5 = rstate.ak4(n);
                const amrex::Real k4 =
                    (rosenbrock::rodas5p_stage_component(state, n) -
                     rstate.y(n) -
                     C::a(5, 1) * k1 -
                     C::a(5, 2) * k2 -
                     C::a(5, 3) * k3) / C::a(5, 4);

                rstate.ak1(n) =
                    C::b(1) * k1 + C::b(2) * k2 + C::b(3) * k3 +
                    C::b(4) * k4 + C::b(5) * k5;
                rstate.ak2(n) =
                    C::c(6, 1) * k1 + C::c(6, 2) * k2 + C::c(6, 3) * k3 +
                    C::c(6, 4) * k4 + C::c(6, 5) * k5;
                rstate.ak3(n) =
                    C::c(7, 1) * k1 + C::c(7, 2) * k2 + C::c(7, 3) * k3 +
                    C::c(7, 4) * k4 + C::c(7, 5) * k5;
                rstate.ak4(n) =
                    C::c(8, 1) * k1 + C::c(8, 2) * k2 + C::c(8, 3) * k3 +
                    C::c(8, 4) * k4 + C::c(8, 5) * k5;
            }

            valid_stage_state = rosenbrock::valid_rodas5p_update_stage(rstate, rstate.ak1);
            if (valid_stage_state) {
                rosenbrock::prepare_rodas5p_update_stage_state(x + C::ct(6) * h,
                                                               state, rstate, rstate.ak1);
                rosenbrock::rodas5p_rhs_from_prepared_state(x + C::ct(6) * h, state, rstate.ak1);
                rstate.n_rhs += 1;
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak1(n) += rstate.ak2(n) / h;
                }
                rosenbrock::solve(rstate, rstate.ak1);
            }

            if (valid_stage_state) {
                for (int n = 1; n <= int_neqs; ++n) {
                    const amrex::Real update =
                        rosenbrock::rodas5p_stage_component(state, n) - rstate.y(n);
                    const amrex::Real k6 = rstate.ak1(n);
                    rstate.ak1(n) = update + k6;
                    rstate.ak3(n) += C::c(7, 6) * k6;
                    rstate.ak4(n) += C::c(8, 6) * k6;
                }
            } else {
                h *= 0.25_rt;
                reject = true;
                n_reject += 1;
                last = false;
                continue;
            }

            valid_stage_state = rosenbrock::valid_rodas5p_update_stage(rstate, rstate.ak1);
            if (valid_stage_state) {
                rosenbrock::prepare_rodas5p_update_stage_state(x + C::ct(7) * h,
                                                               state, rstate, rstate.ak1);
                rosenbrock::rodas5p_rhs_from_prepared_state(x + C::ct(7) * h, state, rstate.ak1);
                rstate.n_rhs += 1;
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak1(n) += rstate.ak3(n) / h;
                }
                rosenbrock::solve(rstate, rstate.ak1);
            }

            if (valid_stage_state) {
                for (int n = 1; n <= int_neqs; ++n) {
                    const amrex::Real update =
                        rosenbrock::rodas5p_stage_component(state, n) - rstate.y(n);
                    const amrex::Real k7 = rstate.ak1(n);
                    rstate.ak1(n) = update + k7;
                    rstate.ak4(n) += C::c(8, 7) * k7;
                }
            } else {
                h *= 0.25_rt;
                reject = true;
                n_reject += 1;
                last = false;
                continue;
            }

            valid_stage_state = rosenbrock::valid_rodas5p_update_stage(rstate, rstate.ak1);
            if (valid_stage_state) {
                rosenbrock::prepare_rodas5p_update_stage_state(x + C::ct(8) * h,
                                                               state, rstate, rstate.ak1);
                rosenbrock::rodas5p_rhs_from_prepared_state(x + C::ct(8) * h, state, rstate.ak1);
                rstate.n_rhs += 1;
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.ak1(n) += rstate.ak4(n) / h;
                }
                rosenbrock::solve(rstate, rstate.ak1);
            }

            if (! valid_stage_state) {
                h *= 0.25_rt;
                reject = true;
                n_reject += 1;
                last = false;
                continue;
            }

            for (int n = 1; n <= int_neqs; ++n) {
                const amrex::Real update =
                    rosenbrock::rodas5p_stage_component(state, n) - rstate.y(n);
                rstate.ak2(n) = rstate.y(n) + update + rstate.ak1(n);
            }

            rstate.n_step += 1;

            const bool valid_state = rosenbrock::valid_integrator_state(rstate.ak2);
            if (! valid_state) {
                reject = true;
                n_reject += 1;
                last = false;
                h *= 0.25_rt;
                continue;
            }

            const bool valid_update = rosenbrock::species_change_valid(rstate, rstate.ak2);
            constexpr amrex::Real err_min = 1.e-10_rt;
            amrex::Real err = rosenbrock::error_norm(rstate, rstate.ak2, rstate.ak1);
            amrex::Real controller_fac = amrex::Clamp(
                std::pow(1.0_rt / amrex::max(err, err_min),
                         1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                std::pow(1.0_rt / amrex::max(errold, err_min),
                         1.0_rt / (integrator_rp::h211b_b * integrator_rp::h211b_k)) *
                std::pow(facold, -1.0_rt / integrator_rp::h211b_b),
                integrator_rp::h211b_fac_min, integrator_rp::h211b_fac_max);
            facold = controller_fac;
            errold = err;

            amrex::Real hnew = h * controller_fac;

            if (err <= 1.0_rt && valid_update) {
                for (int n = 1; n <= int_neqs; ++n) {
                    rstate.y(n) = rstate.ak2(n);
                }
                x += h;

                if (std::abs(hnew) > hmax) {
                    hnew = posneg * hmax;
                }
                if (reject) {
                    hnew = posneg * amrex::min(std::abs(hnew), std::abs(h));
                }
                reject = false;
                n_reject = 0;
                h = hnew;
                break;
            }

            reject = true;
            n_reject += 1;
            last = false;
            if (n_reject >= 2) {
                hnew *= integrator_rp::h211b_reduction_fac;
            }
            if (valid_update) {
                h = posneg * amrex::min(std::abs(hnew),
                                         integrator_rp::h211b_reduction_fac * std::abs(h));
            } else {
                h *= 0.25_rt;
            }
        }
    }
}

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int rosenbrock_integrator (BurnT& state, rosenbrock_t<integrator_neqs<BurnT>()>& rstate)
{
    if (integrator_rp::rosenbrock_tableau == 1) {
        return rosenbrock_integrator<BurnT, rosenbrock::rodas4p_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 2) {
        return rosenbrock_integrator<BurnT, rosenbrock::rodas3p_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 3) {
        return rosenbrock_integrator<BurnT, rosenbrock::ros2s_tableau>(state, rstate);
    }
    if (integrator_rp::rosenbrock_tableau == 4) {
        return rosenbrock_integrator<BurnT, rosenbrock::ros2_tableau>(state, rstate);
    }

    return rosenbrock_integrator<BurnT, rosenbrock::rodas5p_tableau>(state, rstate);
}

#endif
