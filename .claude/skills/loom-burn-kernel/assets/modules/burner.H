// Burner entry point for the Loom burn kernel (see README.md):
// interfaces/burner.H (no NSE), integration/integrator.H,
// integration/VODE/actual_integrator.H, for any rhs.H-framework network + helmholtz + VODE,
// STRANG, HIP branches (see vode.H).
//
// Entry point:
//   MPL_DEVICE void burner(burn_t* state, Real dt, const BurnParams& params, EosTables tab)
// mirrors stock burner(state, dt). BurnParams aggregates the scalar-only
// per-module params structs (read through const references); the helmholtz
// tables stay a by-value EosTables (README rule 8). burn_t is passed by
// pointer because it is written (a mutable struct reference that writes
// scalar fields does not import, probe r15).
#ifndef MPL_BURNER_H
#define MPL_BURNER_H

#include "support.H"
#include "burn.H"
#include "eos.H"
#include "network.H"
#include "vode.H"

namespace mpl {

// All runtime parameters of the burn (scalars only).
struct BurnParams {
    EosParams eos;
    NetParams net;
    IntegratorParams integ;
};

// from integration/VODE/actual_integrator.H:14 actual_integrator
MPL_DEVICE void actual_integrator(burn_t* state, Real dt, bool is_retry,
                                  const IntegratorParams& ip, const EosParams& ep,
                                  const NetParams& np, EosTables tab)
{
    dvode_t vode_state{};
    integrator_setup(state, dt, is_retry, &vode_state, ip, ep, tab);
    state_backup_t state_save = integrator_backup(*state);

    int istate = dvode(state, &vode_state, ip, ep, np, tab);

    integrator_cleanup(&vode_state, state, istate, state_save, dt,
                       vode_final_state_species_failure_tolerance_factor, ip);
}

// from integration/integrator.H:11 integrator_wrapper
template <int enable_retry>
MPL_DEVICE void integrator_wrapper(burn_t* state, Real dt,
                                   const IntegratorParams& ip, const EosParams& ep,
                                   const NetParams& np, EosTables tab)
{
    if constexpr (enable_retry != 0) {
        burn_t old_state = *state;

        actual_integrator(state, dt, false, ip, ep, np, tab);

        if (!state->success) {
            *state = old_state;
            const bool is_retry = true;
            actual_integrator(state, dt, is_retry, ip, ep, np, tab);
        }
    } else {
        actual_integrator(state, dt, false, ip, ep, np, tab);
    }
}

// from integration/integrator.H:33 integrator
MPL_DEVICE void integrator(burn_t* state, Real dt,
                           const IntegratorParams& ip, const EosParams& ep,
                           const NetParams& np, EosTables tab)
{
    if (ip.integrator_use_burn_retry) {
        integrator_wrapper<1>(state, dt, ip, ep, np, tab);
    } else {
        integrator_wrapper<0>(state, dt, ip, ep, np, tab);
    }
}

// from interfaces/burner.H:25 burner (no NSE: integrator(state, dt))
// The params sub-structs are copied to locals once, so everything below
// reads locals rather than the kernel-argument buffer.
MPL_DEVICE void burner(burn_t* state, Real dt, const BurnParams& params, EosTables tab)
{
    const EosParams ep = params.eos;
    const NetParams np = params.net;
    const IntegratorParams ip = params.integ;

    integrator(state, dt, ip, ep, np, tab);
}

} // namespace mpl

#endif
