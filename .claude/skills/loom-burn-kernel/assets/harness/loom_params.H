// Host-side parameter and table setup for the Loom burn kernel (see the
// port's README.md). Shared by the native validators (validate_*.H) and the
// HIP launcher (loom_launch.H). Template: copy to
// unit_test/test_react/loom_params.H and fill in the lines marked "PORT:".
//
// Every fill_* function reads the stock Microphysics globals (runtime
// parameters, helmholtz tables) after eos_init()/network_init() and writes
// the kernel's scalar-only params structs. Bool parameters become int 0/1.
// The EOS ones (fill_eos_params, fill_eos_tables) are network-independent
// and live in <loom_eos_params.H> (copied as is), so the EOS validator builds
// before rates.H/network.H exist. Nothing in either file may be defined again
// in a validator: the validators and the launcher include them.
#ifndef LOOM_PARAMS_H
#define LOOM_PARAMS_H

#include <extern_parameters.H>
#include <AMReX.H>
#include <loom_eos_params.H>  // fill_eos_params, fill_eos_tables

// The kernel directory is LAST on the include path (its headers share names
// with stock ones), so its burner.H is reached by path from this directory
// (LOOM_KERNEL_PATH: <loom_eos_params.H>).
#include LOOM_KERNEL_PATH(burner.H)

// Runtime settings the kernel does not implement: abort with a clear message.
// Called first by fill_burn_params(), so both the launcher and the
// validators reject them.
inline void loom_require_supported_settings()
{
    // Tabulated rates are never ported. PORT: delete this check only if the
    // network has no network.use_tables parameter.
    if (network_rp::use_tables) {
        amrex::Abort("Loom burn kernel: network.use_tables != 0 is not supported "
                     "(tabulated rates are not ported); run with network.use_tables=0");
    }
    // PORT: add a check here for every other runtime setting whose branch
    // the translation omits (abort, never silently ignore).
}

// Rates/screening runtime parameters: one line per field of the port's
// mpl::RatesParams, named <namespace>_<name> after the stock parameter.
inline void fill_rates_params(mpl::RatesParams* p)
{
    // PORT: fill network/rates parameters here, e.g.
    //   p->network_<name> = network_rp::<name> ? 1 : 0;      // bool -> int
    //   p->screening_enable_debye_huckel_skip = screening_rp::enable_debye_huckel_skip ? 1 : 0;
    //   p->screening_debye_huckel_skip_threshold = screening_rp::debye_huckel_skip_threshold;
    (void) p;
}

// Network runtime parameters: every mpl::NetParams field (typically the
// RatesParams ones plus network_use_tables).
inline void fill_net_params(mpl::NetParams* p)
{
    // PORT: fill network parameters here (same lines as fill_rates_params
    // for the shared fields), then:
    p->network_use_tables = network_rp::use_tables ? 1 : 0;
}

// Integrator runtime parameters (+ network_rp::small_x). These are the
// VODE ones; PORT: adjust to the fields of the port's mpl::IntegratorParams
// if it uses another integrator.
inline void fill_integrator_params(mpl::IntegratorParams* p)
{
    p->integrator_rtol_spec = integrator_rp::rtol_spec;
    p->integrator_rtol_enuc = integrator_rp::rtol_enuc;
    p->integrator_atol_spec = integrator_rp::atol_spec;
    p->integrator_atol_enuc = integrator_rp::atol_enuc;
    p->integrator_retry_rtol_spec = integrator_rp::retry_rtol_spec;
    p->integrator_retry_rtol_enuc = integrator_rp::retry_rtol_enuc;
    p->integrator_retry_atol_spec = integrator_rp::retry_atol_spec;
    p->integrator_retry_atol_enuc = integrator_rp::retry_atol_enuc;
    p->integrator_jacobian = integrator_rp::jacobian;
    p->integrator_use_burn_retry = integrator_rp::use_burn_retry ? 1 : 0;
    p->integrator_retry_swap_jacobian = integrator_rp::retry_swap_jacobian ? 1 : 0;
    p->integrator_scale_system = integrator_rp::scale_system ? 1 : 0;
    p->integrator_call_eos_in_rhs = integrator_rp::call_eos_in_rhs ? 1 : 0;
    p->integrator_integrate_energy = integrator_rp::integrate_energy ? 1 : 0;
    p->integrator_use_number_densities = integrator_rp::use_number_densities ? 1 : 0;
    p->integrator_subtract_internal_energy = integrator_rp::subtract_internal_energy ? 1 : 0;
    p->integrator_do_species_clip = integrator_rp::do_species_clip ? 1 : 0;
    p->integrator_renormalize_abundances = integrator_rp::renormalize_abundances ? 1 : 0;
    p->integrator_do_corrector_validation = integrator_rp::do_corrector_validation ? 1 : 0;
    p->integrator_correct_jacobian_for_const_e = integrator_rp::correct_jacobian_for_const_e ? 1 : 0;
    p->integrator_ode_max_steps = integrator_rp::ode_max_steps;
    p->integrator_ode_max_dt = integrator_rp::ode_max_dt;
    p->integrator_SMALL_X_SAFE = integrator_rp::SMALL_X_SAFE;
    p->integrator_MAX_TEMP = integrator_rp::MAX_TEMP;
    p->integrator_react_boost = integrator_rp::react_boost;
    p->integrator_species_failure_tolerance = integrator_rp::species_failure_tolerance;
    p->integrator_X_reject_buffer = integrator_rp::X_reject_buffer;
    p->network_small_x = network_rp::small_x;
}

// Everything the burn reads. Aborts on unsupported settings.
inline void fill_burn_params(mpl::BurnParams* p)
{
    loom_require_supported_settings();
    fill_eos_params(&p->eos);
    fill_net_params(&p->net);
    fill_integrator_params(&p->integ);
    // PORT: fill any other mpl::BurnParams members the port's burner adds.
}

#endif
