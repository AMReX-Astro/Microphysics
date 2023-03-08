#ifndef INTEGRATOR_RHS_STRANG_H
#define INTEGRATOR_RHS_STRANG_H

#include <network.H>
#include <actual_network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_rhs.H>
#endif
#include <burn_type.H>
#include <extern_parameters.H>
#include <integrator_data.H>
#include <integrator_type_strang.H>
#ifdef NONAKA_PLOT
#include <nonaka_plot.H>
#endif

// The rhs routine provides the right-hand-side for the DVODE solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs (const Real time, burn_t& state, T& int_state, RArray1D& ydot, [[maybe_unused]] const bool in_jacobian=false)
{

    // We are integrating a system of
    //
    // y(1:NumSpec) = dX/dt
    // y(net_ienuc) = denuc/dt

    // we come in with (integrator type) int_state having the current
    // solution (or predicted) state and we need to copy this over to
    // the (burn_t) state to interface with the actual reaction
    // network

    // Fix the state as necessary -- this ensures that the mass
    // fractions that enter are valud (and optionally normalized)

    clean_state(int_state);

    // Update the thermodynamics as necessary -- this primarily takes
    // the information in int_state (X and e), copies it to the
    // (burn_t) state and calls the EOS to get state.T

    update_thermodynamics(state, int_state);

    // Only do the burn if the incoming temperature is within the temperature
    // bounds. Otherwise set the RHS to zero and return.

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int n = 1; n <= INT_NEQS; ++n) {
            ydot(n) = 0.0_rt;
        }

        return;

    }

    state.time = time;

    // at this point, the (burn_t) state is synchronized with the integrator
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
    // turn it off for primordial chem
    if (!use_number_densities) {
        for (int n = 1; n <= NumSpec; ++n) {
            ydot(n) *= aion[n-1];
        }
    }

    // Allow energy integration to be disabled.

    if (!integrate_energy) {
        ydot(net_ienuc) = 0.0_rt;
    }

    // apply fudge factor:

    if (react_boost > 0.0_rt) {
        for (int n = 1; n <= INT_NEQS; ++n) {
            ydot(n) *= react_boost;
        }
    }

}



// Analytical Jacobian
template<typename T, class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac ([[maybe_unused]] const Real time, burn_t& state, T& int_state, MatrixType& pd)
{

    // Only do the burn if the incoming temperature is within the temperature
    // bounds. Otherwise set the Jacobian to zero and return.

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int j = 1; j <= INT_NEQS; ++j) {
            for (int i = 1; i <= INT_NEQS; ++i) {
                pd(i,j) = 0.0_rt;
            }
        }

        return;

    }

    // Call the specific network routine to get the Jacobian.

    integrator_to_burn(int_state, state);

#ifdef NEW_NETWORK_IMPLEMENTATION
    RHS::jac(state, pd);
#else
    actual_jac(state, pd);
#endif

    // We integrate X, not Y
    // turn it off for primordial chem
    if (!use_number_densities) {
        for (int j = 1; j <= NumSpec; ++j) {
            for (int i = 1; i <= INT_NEQS; ++i) {
                pd.mul(j, i, aion[j-1]);
                pd.mul(i, j, aion_inv[j-1]);
            }
        }
    }

    // apply fudge factor:
    if (react_boost > 0.0_rt) {
        pd.mul(react_boost);
    }

    // Allow temperature and energy integration to be disabled.

    if (!integrate_energy) {
        for (int j = 1; j <= INT_NEQS; ++j) {
            pd(net_ienuc,j) = 0.0_rt;
        }
    }

}

#endif
