#ifndef RHS_ZONES_H
#define RHS_ZONES_H

#include <variables.H>
#include <network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#endif
#include <eos_type.H>
#include <eos.H>
#include <burn_type.H>
#include <burner.H>
#include <extern_parameters.H>
#include <ArrayUtilities.H>
#include <variables.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool do_rhs (int i, int j, int k, Array4<Real> const& state, const plot_t& p)
{

    burn_t burn_state;

    burn_state.rho = state(i, j, k, p.irho);
    burn_state.T = state(i, j, k, p.itemp);
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = state(i, j, k, p.ispec_old + n);
    }

    normalize_abundances_burn(burn_state);

    eos_t eos_state;

    burn_to_eos(burn_state, eos_state);

    eos(eos_input_rt, eos_state);

    eos_to_burn(eos_state, burn_state);

    // the integrator doesn't actually care about the initial internal
    // energy.
    burn_state.e = 0.0_rt;

    Array1D<Real, 1, neqs> ydot;

    ArrayUtil::MathArray2D<1, neqs, 1, neqs> jac;

#ifdef NEW_NETWORK_IMPLEMENTATION
    RHS::rhs(burn_state, ydot);
    RHS::jac(burn_state, jac);
#else
    actual_rhs(burn_state, ydot);
    actual_jac(burn_state, jac);
#endif

    for (int n = 1; n <= NumSpec; ++n) {
        state(i, j, k, p.ispec + n - 1) = ydot(n);
        state(i, j, k, p.ienuc_dot) = ydot(net_ienuc);
    }

    int n = 0;
    for (int jj = 1; jj <= neqs; ++jj) {
        for (int ii = 1; ii <= neqs; ++ii) {
            state(i, j, k, p.ijac + n) = jac.get(ii, jj);
            n++;
        }
    }

    return true;

}

#endif
