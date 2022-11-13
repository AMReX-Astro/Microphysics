#ifndef _vode_rhs_H_
#define _vode_rhs_H_

#include <network.H>
#include <actual_network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_rhs.H>
#endif
#include <burn_type.H>
#include <extern_parameters.H>
#include <vode_type.H>
#ifdef NONAKA_PLOT
#include <nonaka_plot.H>
#endif

// The rhs routine provides the right-hand-side for the DVODE solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs (const Real time, burn_t& state, dvode_t& vode_state, RArray1D& ydot, const bool in_jacobian=false)
{

    // We are integrating a system of
    //
    // y(1:NumSpec) = dX/dt
    // y(net_ienuc) = denuc/dt

    // we come in with (dvode_t) vode_state having the current
    // solution (or predicted) state and we need to copy this over to
    // the (burn_t) state to interface with the actual reaction
    // network

    // Fix the state as necessary -- this ensures that the mass
    // fractions that enter are valud (and optionally normalized)

    clean_state(vode_state);

    // Update the thermodynamics as necessary -- this primarily takes
    // the information in vode_state (X and e), copies it to the
    // (burn_t) state and calls the EOS to get state.T

    update_thermodynamics(state, vode_state);

    // Only do the burn if the incoming temperature is within the temperature
    // bounds. Otherwise set the RHS to zero and return.

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int n = 1; n <= VODE_NEQS; ++n) {
            ydot(n) = 0.0_rt;
        }

        return;

    }

    state.time = vode_state.t;

    // at this point, the (burn_t) state is synchronized with VODE
    // and is thermodynamically consistent.

    // Call the specific network routine to get the RHS.

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

    // We integrate X, not Y

    for (int n = 1; n <= NumSpec; ++n) {
        ydot(n) *= aion[n-1];
    }

    // Allow energy integration to be disabled.

    if (!integrate_energy) {
        ydot(net_ienuc) = 0.0_rt;
    }

    // apply fudge factor:

    if (react_boost > 0.0_rt) {
        for (int n = 1; n <= VODE_NEQS; ++n) {
            ydot(n) *= react_boost;
        }
    }

}



// Analytical Jacobian
template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (burn_t& state, dvode_t& vode_state, MatrixType& pd)
{
    // NOTE: the time at which to evaluate the Jacobian is not
    // explicitly passed. VODE always evaluates the analytic
    // Jacobian at vode_state.tn, and we pass vode_state.

    // Only do the burn if the incoming temperature is within the temperature
    // bounds. Otherwise set the Jacobian to zero and return.

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int j = 1; j <= VODE_NEQS; ++j) {
            for (int i = 1; i <= VODE_NEQS; ++i) {
                pd(i,j) = 0.0_rt;
            }
        }

        return;

    }

    // Call the specific network routine to get the Jacobian.

    vode_to_burn(vode_state, state);

#ifdef NEW_NETWORK_IMPLEMENTATION
    RHS::jac(state, pd);
#else
    actual_jac(state, pd);
#endif

    // We integrate X, not Y
    for (int j = 1; j <= NumSpec; ++j) {
        for (int i = 1; i <= VODE_NEQS; ++i) {
            pd.mul(j, i, aion[j-1]);
            pd.mul(i, j, aion_inv[j-1]);
        }
    }

    // apply fudge factor:
    if (react_boost > 0.0_rt) {
        pd.mul(react_boost);
    }

    // Allow temperature and energy integration to be disabled.

    if (!integrate_energy) {
        for (int j = 1; j <= VODE_NEQS; ++j) {
            pd(net_ienuc,j) = 0.0_rt;
        }
    }

}

#endif
