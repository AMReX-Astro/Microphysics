#ifndef VODE_UTIL_H
#define VODE_UTIL_H

// Common variables and routines for burners
// that use VODE for their integration.

#include <extern_parameters.H>
#include <vode_type.H>
#include <vode_dvode.H>

template <typename IntT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void setup_integrator_type (IntT& vode_state) {

    // set the Jacobian type

    vode_state.jacobian_type = jacobian;

    // Set the (inverse of the) timestep limiter.

    vode_state.HMXI = 1.0_rt / integrator_rp::ode_max_dt;

}

template <typename BurnT, typename IntT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void teardown_integrator_type (BurnT& state, const IntT& vode_state) {

    // Get the number of RHS and Jacobian evaluations.

    state.n_rhs = vode_state.NFE;
    state.n_jac = vode_state.NJE;
    state.n_step = vode_state.NST;

}

template <typename BurnT, typename IntT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int do_integration (BurnT& state, IntT& vode_state) {

    return dvode(state, vode_state);
}

#endif
