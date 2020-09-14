#ifndef VODE_RHS_SIMPLIFIED_SDC_H
#define VODE_RHS_SIMPLIFIED_SDC_H

#include <network.H>
#include <actual_network.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <extern_parameters.H>
#include <vode_type.H>

// The f_rhs routine provides the right-hand-side for the DVODE solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void f_rhs(const Real time, burn_t& state, dvode_t& vode_state, RArray1D& ydot)
{

  // Fix the state as necessary.

  clean_state(time, vode_state);

  // convert to the burn_t

  vode_to_burn(time, vode_state, state);

  state.time = time;

  // call the specific network to get the RHS
  YdotNetArray1D ydot_react = {0};

  actual_rhs(state, ydot_react);

  // apply fudge factor:
  if (react_boost > 0.0_rt) {
    for (int n = 1; n <= VODE_NEQS; ++n) {
      ydot_react(n) *= react_boost;
    }
  }

  // convert back to the vode type -- this will add the advective terms

  rhs_to_vode(time, state, ydot_react, vode_state, ydot);

}



// Analytical Jacobian

template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (burn_t& state, dvode_t& vode_state, MatrixType& pd)
{
    // NOTE: the time at which to evaluate the Jacobian is not
    // explicitly passed. VODE always evaluates the analytic
    // Jacobian at vode_state.tn, and we pass vode_state.

    // Call the specific network routine to get the Jacobian.

    vode_to_burn(vode_state.tn, vode_state, state);

    JacNetArray2D jac_react;

    actual_jac(state, jac_react);

    // apply fudge factor:

    if (react_boost > 0.0_rt) {
        for (int j = 1; j <= VODE_NEQS; ++j) {
            for (int i = 1; i <= VODE_NEQS; ++i) {
                jac_react.mul(i, j, react_boost);
            }
        }
    }

    // convert to the system we are using

    jac_to_vode(vode_state.tn, jac_react, vode_state, pd);

}
#endif
