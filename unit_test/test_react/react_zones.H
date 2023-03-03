#ifndef REACT_ZONES_H
#define REACT_ZONES_H

#include <variables.H>
#include <network.H>
#include <eos_type.H>
#include <eos.H>
#include <burn_type.H>
#include <burner.H>
#include <extern_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool do_react (int i, int j, int k, Array4<Real> const& state,
               burn_t& burn_state, Array4<int> const& n_rhs, const plot_t& p)
{

    burn_state.rho = state(i, j, k, p.irho);
    burn_state.T = state(i, j, k, p.itemp);
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = state(i, j, k, p.ispec_old + n);
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
      burn_state.aux[n] = state(i, j, k, p.iaux_old + n);
    }
#endif

    normalize_abundances_burn(burn_state);

    // the integrator doesn't actually care about the initial internal
    // energy.
    burn_state.e = 0.0_rt;

    Real dt = tmax;

    burn_state.i = i;
    burn_state.j = j;
    burn_state.k = k;
#ifdef NONAKA_PLOT
    burn_state.level = 0;
    burn_state.reference_time = 0.0;
#endif

    burn_state.T_fixed = -1.0_rt;

    burner(burn_state, dt);

    for (int n = 0; n < NumSpec; ++n) {
        state(i, j, k, p.ispec + n) = burn_state.xn[n];
        state(i, j, k, p.irodot + n) = (burn_state.xn[n] - state(i, j, k, p.ispec_old + n)) / dt;
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
      state(i, j, k, p.iaux + n) = burn_state.aux[n];
    }
#endif

    state(i, j, k, p.irho_hnuc) = state(i, j, k, p.irho) * burn_state.e / dt;

    n_rhs(i, j, k, 0) = burn_state.n_rhs;
    n_rhs(i, j, k, 1) = burn_state.n_step;

    return burn_state.success;

}

#endif
