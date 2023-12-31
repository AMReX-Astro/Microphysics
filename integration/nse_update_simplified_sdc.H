#ifndef NSE_UPDATE_H
#define NSE_UPDATE_H

#include <iostream>
#include <fstream>
#include <actual_network.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>

#include <burn_type.H>
#include <nse_eos.H>

#ifdef NSE_TABLE
#include <nse_table.H>
#endif
#ifdef NSE_NET
#include <nse_solver.H>
#endif

using namespace amrex;

#if defined(NSE_TABLE)

///
/// this acts as an explicit Euler step for the system (rho e, rho aux)
/// on input, *_source are the reactive sources at time t0 and on output
/// they are the sources at time t0+dt
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_derivs(const Real rho0, const Real rhoe0, const Real *rhoaux0,
                const Real dt, const Real *ydot_a,
                Real& drhoedt, Real* drhoauxdt, const Real T_fixed) {

    // start with the current state and compute T

    Real T0;
    Real abar;
    Real Ye0 = rhoaux0[iye] / rho0;

    if (T_fixed > 0) {
        T0 = T_fixed;
        abar = rhoaux0[iabar] / rho0;
    } else {
        Real e0 = rhoe0 / rho0;
        T0 = 1.e8;  // initial guess
        nse_T_abar_from_e(rho0, e0, Ye0, T0, abar);
    }

    // compute the plasma neutrino losses at t0

    Real snu{0.0};
    Real dsnudt{0.0};
    Real dsnudd{0.0};
    Real dsnuda{0.0};
    Real dsnudz{0.0};

    Real zbar = abar * Ye0;

#ifdef NEUTRINOS
    constexpr int do_derivatives = 0;
    sneut5<do_derivatives>(T0, rho0, abar, zbar,
                           snu, dsnudt, dsnudd, dsnuda, dsnudz);
#endif

    // call the NSE table at t0

    constexpr bool skip_X_fill{true};

    nse_table_t nse_state;
    nse_state.T = T0;
    nse_state.rho = rho0;
    nse_state.Ye = Ye0;
    nse_interp(nse_state, skip_X_fill);

    Real abar0_out = nse_state.abar;
    Real bea0_out = nse_state.bea;
    Real dyedt0 = nse_state.dyedt;

    // construct initial sources at t0

    Real rhoe_source{};
    if (integrator_rp::nse_include_enu_weak == 1) {
        rhoe_source = -rho0 * (nse_state.e_nu + snu);
    } else {
        rhoe_source = -rho0 * snu;
    }

    Real rhoaux_source[NumAux];
    rhoaux_source[iye] = rho0 * dyedt0;
    rhoaux_source[iabar] = 0.0;
    rhoaux_source[ibea] = rho0 * nse_state.dbeadt;

    // evolve for eps * dt;

    Real tau = integrator_rp::nse_deriv_dt_factor * dt;

    Real rho1 = rho0 + tau * ydot_a[SRHO];
    Real rhoe1 = rhoe0 + tau * ydot_a[SEINT] + tau * rhoe_source;
    Real rhoaux1[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        rhoaux1[n] = rhoaux0[n] + tau * ydot_a[SFX+n] + tau * rhoaux_source[n];
    }

    // compute the temperature at t0 + tau

    Real T1;
    Real abar1_out;
    Real Ye1 = rhoaux1[iye] / rho1;

    if (T_fixed > 0) {
        T1 = T_fixed;
        abar1_out = rhoaux1[iabar] / rho1;
    } else {
        Real e1 = rhoe1 / rho1;
        T1 = T0;
        nse_T_abar_from_e(rho1, e1, Ye1, T1, abar1_out);
    }

    // call NSE at t0 + tau

    nse_state.T = T1;
    nse_state.rho = rho1;
    nse_state.Ye = Ye1;

    nse_interp(nse_state, skip_X_fill);

    Real bea1_out = nse_state.bea;

    // update abar -- this will be the same as nse_T_abar_from_e, but
    // for the case with T_fixed > 0, this abar will be consistent
    // with NSE.
    abar1_out = nse_state.abar;

    // construct the finite-difference approximation to the derivatives

    // note that abar, (B/A) here have already seen advection
    // implicitly

    Real rho_bea_tilde = rho1 * bea1_out - tau * ydot_a[SFX+ibea];
    Real rho_dBEA = rho_bea_tilde - rho0 * bea0_out; // this is MeV / nucleon * g / cm**3

    Real rho_abar_tilde = rho1 * abar1_out - tau * ydot_a[SFX+iabar];
    Real rho_dabar = rho_abar_tilde - rho0 * abar0_out; // this is MeV / nucleon * g / cm**3

    drhoedt = rho_dBEA * C::MeV2eV * C::ev2erg * C::n_A / tau;
    if (integrator_rp::nse_include_enu_weak == 1) {
        drhoedt -= rho0 * (nse_state.e_nu + snu);
    } else {
        drhoedt -= rho0 * snu;
    }
    drhoauxdt[iabar] = rho_dabar / tau;
    drhoauxdt[iye] = rho0 * dyedt0;
    drhoauxdt[ibea] = rho_dBEA / tau;

}

///
/// update the state due to NSE changes for the simplified-SDC case
/// this version works with the tabulated NSE and requires AUX_THERMO
///
template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sdc_nse_burn(BurnT& state, const Real dt) {

    using namespace AuxZero;

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // store the initial state

    Real rho_old = state.y[SRHO];
    Real rhoe_old = state.y[SEINT];
    Real rhoaux_old[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        rhoaux_old[n] = state.y[SFX+n];
    }

    // density and momentum have no reactive sources

    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    // do an RK2 integration

    // get the derivatives at t = t^n

    Real drhoedt;
    Real drhoauxdt[NumSpec];

    nse_derivs(rho_old, rhoe_old, rhoaux_old,
               dt, state.ydot_a,
               drhoedt, drhoauxdt, state.T_fixed);

    // evolve to the midpoint in time

    Real rho_tmp = rho_old + 0.5_rt * dt * state.ydot_a[SRHO];
    Real rhoe_tmp = rhoe_old + 0.5_rt * dt * (state.ydot_a[SEINT] + drhoedt);
    Real rhoaux_tmp[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        rhoaux_tmp[n] = rhoaux_old[n] +
            0.5_rt * dt * (state.ydot_a[SFX+n] + drhoauxdt[n]);
    }

    // compute the derivatives at the midpoint in time

    nse_derivs(rho_tmp, rhoe_tmp, rhoaux_tmp,
               dt, state.ydot_a,
               drhoedt, drhoauxdt, state.T_fixed);

    // evolve to the new time

    Real rho_new = rho_old + dt * state.ydot_a[SRHO];
    Real rhoe_new = rhoe_old + dt * (state.ydot_a[SEINT] + drhoedt);
    Real rhoaux_new[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        rhoaux_new[n] = rhoaux_old[n] +
            dt * (state.ydot_a[SFX+n] + drhoauxdt[n]);
    }

    // get the new temperature.  This will also get a new abar, but
    // we'll use the version from the NSE table compute below for the
    // final abar.

    Real T_new;
    Real abar_new;
    Real Ye_new = rhoaux_new[iye] / rho_new;

    if (state.T_fixed > 0) {
        T_new = state.T_fixed;
        // note: this is not consistent with NSE
        abar_new = rhoaux_new[iabar] / rho_new;
    } else {
        Real e_new = rhoe_new / rho_new;
        T_new = 1.e8;  // initial guess
        nse_T_abar_from_e(rho_new, e_new, Ye_new, T_new, abar_new);
    }

    // do a final NSE call -- we want the ending B/A to be consistent
    // with our thermodynamic state here, plus we need to get the mass
    // fractions.  We'll also use the abar from this call.  This should
    // be the same as what we got from nse_T_abar_from_e, but if we are
    // doing T_fixed, this abar will be consistent with NSE.


    constexpr bool skip_X_fill{false};

    nse_table_t nse_state;
    nse_state.T = T_new;
    nse_state.rho = rho_new;
    nse_state.Ye = Ye_new;

    nse_interp(nse_state, skip_X_fill);

    // store the new state

    // the new mass fractions are just those that come from the table
    // make sure they are normalized

    Real sum_X{0.0_rt};
    for (auto & xn : nse_state.X) {
        xn = std::clamp(xn, small_x, 1.0_rt);
        sum_X += xn;
    }

    for (auto & xn : nse_state.X) {
        xn /= sum_X;
    }

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = rho_new * nse_state.X[n];
    }

    // aux data comes from the integration or the table

    state.y[SFX+iye] = rhoaux_new[iye];
    state.y[SFX+iabar] = rho_new * nse_state.abar;
    state.y[SFX+ibea] = rho_new * nse_state.bea;

    // density and momenta have already been updated

    // get the energy release from the change in rhoe over the timestep
    // excluding any advection, and use that to update the total energy

    Real rho_enucdot = (rhoe_new - rhoe_old) / dt - state.ydot_a[SEINT];

    state.y[SEINT] = rhoe_new;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_enucdot;

}

#else
///
/// update the state due to NSE changes for the simplified-SDC case
/// this version works with the self-consistent NSE
///
template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sdc_nse_burn(BurnT& state, const Real dt) {

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // store the initial mass fractions -- we will need these
    // to compute the energy release.

    Real X_old[NumSpec];

    for (int n = 0; n < NumSpec; ++n) {
        X_old[n] = state.y[SFS+n] / state.y[SRHO];
    }

    // density and momentum have no reactive sources
    Real rho_old = state.y[SRHO];

    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    // get the neutrino loss term -- we want to use the state that we
    // came in here with, so the original Abar and Zbar
    Real snu{0.0};
    Real dsnudt{0.0};
    Real dsnudd{0.0};
    Real dsnuda{0.0};
    Real dsnudz{0.0};

    Real abar{0.0};
    Real zbar{0.0};
    for (int n = 0; n < NumSpec; ++n) {
        abar += X_old[n] * aion_inv[n];
        zbar += X_old[n] * zion[n] * aion_inv[n];
    }
    abar = 1.0 / abar;
    zbar *= abar;

#ifdef NEUTRINOS
    constexpr int do_derivatives = 0;
    sneut5<do_derivatives>(T_in, rho_old, abar, zbar,
                           snu, dsnudt, dsnudd, dsnuda, dsnudz);
#endif
    Real snu_old = snu;


    // if our network could return the evolution of Ye due to the
    // weak interactions, we would evaluate the NSE state here and
    // get dYe/dt.

    Real dyedt_old = 0.0;


    // predict the U^{n+1,*} state with only estimates of the X
    // updates with advection to dt and the neutrino loss term in
    // energy

    BurnT burn_state;
    burn_state.T = T_in; // initial guess
    burn_state.mu_p = state.mu_p;
    burn_state.mu_n = state.mu_n;

    Real rhoX_source[NumSpec] = {0.0_rt};

    Real rhoX_new[NumSpec];

    Real rhoe_new;
    Real rho_enucdot = -rho_old * snu;

    Real rho_half = 0.5_rt * (rho_old + state.y[SRHO]);

    // predict the new aux for the first iteration -- this is really
    // just including the advection bits

    for (int n = 0; n < NumSpec; n++) {
        rhoX_new[n] = state.y[SFS+n] + dt * state.ydot_a[SFS+n] + dt * rhoX_source[n];
    }

    burn_t nse_state;

    for (int iter = 0; iter < integrator_rp::nse_iters; iter++) {

        // update (rho e)^{n+1} based on the new energy generation rate
        rhoe_new = state.y[SEINT] + dt * state.ydot_a[SEINT] + dt * rho_enucdot;

        // call the EOS to get the updated T*

        Real T_new;
        burn_state.rho = state.y[SRHO];
        burn_state.e = rhoe_new / state.y[SRHO];
        for (int n = 0; n < NumSpec; n++) {
            burn_state.xn[n] = rhoX_new[n] / state.y[SRHO];
            burn_state.y[SFS+n] = rhoX_new[n];
        }

        if (state.T_fixed > 0) {
            T_new = state.T_fixed;
        } else {
            eos(eos_input_re, burn_state);
            T_new = burn_state.T;
        }
        burn_state.T = T_new;

        // solve for the NSE state for this network we need to call
        // get_nse_state with a burn_t.  We will have it compute
        // Ye from the input X's

        nse_state = get_actual_nse_state(burn_state);

        // compute the energy release.  The mass fractions in nse_state.xn[]
        // include the advective parts, so first we need to remove that.

        Real rhoX_tilde[NumSpec];
        for (int n = 0; n < NumSpec; ++n) {
            rhoX_tilde[n] = state.y[SRHO] * nse_state.xn[n] - dt * state.ydot_a[SFS+n];
        }

        Real dyedt = 0.0_rt;   // we can update this in the future by calling actual_rhs()

        // we want to compute (rho eps) = - N_A c^2 sum{m_i (rhoX_tilde - rhoX_old) / A_i}
        rho_enucdot = 0.0;
        for (int n = 0; n < NumSpec; ++n) {
            rho_enucdot += (rhoX_tilde[n] - rho_old * X_old[n]) *
                network::mion(n+1) * aion_inv[n];
        }
        rho_enucdot *= C::Legacy::enuc_conv2;

        // now get the updated neutrino term
        abar = 0.0;
        zbar = 0.0;
        for (int n = 0; n < NumSpec; ++n) {
            abar += nse_state.xn[n] * aion_inv[n];
            zbar += nse_state.xn[n] * zion[n] * aion_inv[n];
        }
        abar = 1.0 / abar;
        zbar *= abar;

#ifdef NEUTRINOS
        constexpr int do_derivatives = 0;
        sneut5<do_derivatives>(T_new, state.y[SRHO], abar, zbar,
                               snu, dsnudt, dsnudd, dsnuda, dsnudz);
#endif

        rho_enucdot -= 0.5_rt * rho_half * (snu_old + snu);

        // update the new state for the next pass -- this should
        // already implicitly have the advective portion included,
        // since it was there when we called the NSE state

        for (int n = 0; n < NumSpec; n++) {
            rhoX_new[n] = state.y[SRHO] * nse_state.xn[n];
        }

    }

    // now update the aux quantities

    // the new mass fractions are just those that come from the table
    // make sure they are normalized
    Real sum_X{0.0_rt};
    Real X[NumSpec] = {0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        X[n] = amrex::max(small_x, amrex::min(1.0_rt, nse_state.xn[n]));
        sum_X += X[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        X[n] /= sum_X;
    }

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = state.y[SRHO] * X[n];
    }


    // density and momenta have already been updated

    // update the total and internal energy now

    state.y[SEINT] += dt * state.ydot_a[SEINT] + dt * rho_enucdot;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_enucdot;

    // store the chemical potentials
    state.mu_p = nse_state.mu_p;
    state.mu_n = nse_state.mu_n;

}
#endif

#endif
