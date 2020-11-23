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
bool do_react (const plot_t p, const int i, const int j, const int k,
               Array4<Real> const& state, Array4<int> const& n_rhs)
{

    burn_t burn_state;

    // populate the SDC state.  Our strategy for the unit test is to choose
    // the advective terms to all be zero and to choose the velocity to be
    // Mach = 0.1

    // call the EOS first to get the sound speed

    eos_t eos_state;
    eos_state.rho = state(i, j, k, p.irho);
    eos_state.T = state(i, j, k, p.itemp);
    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = state(i, j, k, p.ispec_old+n);
    }
#ifdef NSE_THERMO
    set_nse_aux_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

    // now fill the burn state

    burn_state.rho = eos_state.rho;
    burn_state.T = eos_state.T;

    for (int n = 0; n < NumSpec; n++) {
      burn_state.y[SFS+n] = burn_state.rho * eos_state.xn[n];
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        burn_state.y[SFX+n] = burn_state.rho * eos_state.aux[n];
    }
#endif

#if defined(SDC_EVOLVE_ENERGY)
    burn_state.y[SRHO] = state(i, j, k, p.irho);

    // we will pick velocities to be 10% of the sound speed
    burn_state.y[SMX] = burn_state.y[SRHO] * 0.1_rt * eos_state.cs;
    burn_state.y[SMY] = burn_state.y[SRHO] * 0.1_rt * eos_state.cs;
    burn_state.y[SMZ] = burn_state.y[SRHO] * 0.1_rt * eos_state.cs;

    burn_state.y[SEINT] = burn_state.y[SRHO] * eos_state.e;
    burn_state.y[SEDEN] = burn_state.y[SEINT] +
      0.5_rt * (burn_state.y[SMX] * burn_state.y[SMX] +
                burn_state.y[SMY] * burn_state.y[SMY] +
                burn_state.y[SMZ] * burn_state.y[SMZ]) / burn_state.y[SRHO];

    // get T from total energy
    burn_state.T_from_eden = true;

#elif defined(SDC_EVOLVE_ENTHALPY)
    burn_state.y[SENTH] = burn_state.rho * eos_state.h;

    burn_state.p0 = eos_state.p;
#endif

    // normalize
    Real sum_spec = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
      sum_spec += burn_state.y[SFS+n] / burn_state.rho;
    }
    for (int n = 0; n < NumSpec; n++) {
      burn_state.y[SFS+n] /= sum_spec;
    }

#if defined(SDC_EVOLVE_ENERGY)
    Real rhoe_old = burn_state.y[SEINT];
#elif defined (SDC_EVOLVE_ENTHALPY)
    Real rhoh_old = burn_state.y[SENTH];
#endif

    // zero out the advective terms

    for (int n = 0; n < SVAR; n++) {
        burn_state.ydot_a[n] = 0.0;
    }

    // these need to be initialized 

    burn_state.sdc_iter = 1;
    burn_state.num_sdc_iters = 1;

    integrator(burn_state, tmax);

    for (int n = 0; n < NumSpec; ++n) {
        state(i, j, k, p.ispec + n) = burn_state.y[SFS+n]/burn_state.rho;
        state(i, j, k, p.irodot + n) = (burn_state.y[SFS+n] -
                                        state(i, j, k, p.irho) * state(i, j, k, p.ispec_old + n)) / tmax;
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        state(i, j, k, p.iaux + n) = burn_state.y[SFX+n]/burn_state.rho;
    }
#endif

#if defined(SDC_EVOLVE_ENERGY)
    state(i, j, k, p.irho_Hnuc) = (burn_state.y[SEINT] - rhoe_old) / tmax;
#elif defined (SDC_EVOLVE_ENTHALPY)
    state(i, j, k, p.irho_Hnuc) = (burn_state.y[SENTH] - rhoh_old) / tmax;
#endif

    n_rhs(i, j, k, 0) = burn_state.n_rhs;

    return burn_state.success;

}

#endif
