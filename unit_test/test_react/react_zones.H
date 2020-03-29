#ifndef _react_zones_H_
#define _react_zones_H_

#include <variables.H>
#include <network.H>
#include <eos_type.H>
#include <eos.H>
#include <burn_type.H>
#include <burner.H>
#include <extern_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void do_react (int i, int j, int k, Array4<Real> const& state, Array4<int> const& n_rhs)
{

    burn_t burn_state;

    burn_state.rho = state(i, j, k, p_irho);
    burn_state.T = state(i, j, k, p_itemp);
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = state(i, j, k, p_ispec_old + n);
    }

    normalize_abundances_burn(burn_state);

    // the integrator doesn't actually care about the initial internal
    // energy.
    burn_state.e = 0.0_rt;

    Real dt = tmax;

    burner(burn_state, dt);

    for (int n = 0; n < NumSpec; ++n) {
        state(i, j, k, p_ispec + n) = burn_state.xn[n];
        state(i, j, k, p_irodot + n) = (burn_state.xn[n] - state(i, j, k, p_ispec_old + n)) / dt;
    }

    state(i, j, k, p_irho_hnuc) = state(i, j, k, p_irho) * burn_state.e / dt;

    n_rhs(i, j, k, 0) = burn_state.n_rhs;

}

#endif
