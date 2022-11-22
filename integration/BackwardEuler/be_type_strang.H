#ifndef BE_TYPE_STRANG_H
#define BE_TYPE_STRANG_H

#include <be_type.H>
#include <network.H>
#include <actual_network.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>
#include <be_integrator.H>

#define VODELIKE_ERROR 1

using namespace integrator_rp;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state (burn_t& state)
{
    // Renormalize the abundances.

    normalize_abundances_burn(state);

    // Evaluate the EOS to get T from e.

    if (call_eos_in_rhs) {
        eos(eos_input_re, state);
    }

    // Ensure that the temperature always stays within reasonable limits.

    state.T = amrex::min(MAX_TEMP, amrex::max(state.T, EOSData::mintemp));

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void evaluate_rhs(burn_t& state, const be_t& be, Array1D<Real, 1, BE_NEQS>& ydot)
{

    for (int ii = 1; ii <= NumSpec; ii++) {
        state.xn[ii-1] = be.y(ii);
    }
    state.e = be.y(net_ienuc);

    clean_state(state);

    actual_rhs(state, ydot);

    // convert species to be in terms of mass fractions

    for (int ii = 1; ii <= NumSpec; ++ii) {
        ydot(ii) *= aion[ii-1];
    }

    // if we are not integrating e, then zero this out

    if (!integrate_energy) {
        ydot(net_ienuc) = 0.0_rt;
    }
}


#endif
