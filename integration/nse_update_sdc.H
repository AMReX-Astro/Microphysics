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

#ifdef NSE_TABLE
#include <nse_table.H>
#include <nse_eos.H>
#endif
#ifdef NSE_NET
#include <nse_solver.H>
#include <nse_eos.H>
#endif

using namespace amrex::literals;

#if defined(NSE_TABLE)

///
/// this acts as an explicit Euler step for the system (rho e, rho aux)
/// on input, *_source are the reactive sources at time t0 and on output
/// they are the sources at time t0+dt
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_derivs(const amrex::Real rho0, const amrex::Real rhoe0, const amrex::Real *rhoaux0,
                const amrex::Real dt, const amrex::Real *ydot_a,
                amrex::Real& drhoedt, amrex::Real* drhoauxdt, const amrex::Real T_fixed) {

    // start with the current state and compute T

    amrex::Real T0;
    amrex::Real abar;
    amrex::Real Ye0 = rhoaux0[iye] / rho0;

    if (T_fixed > 0) {
        T0 = T_fixed;
        abar = rhoaux0[iabar] / rho0;
    } else {
        amrex::Real e0 = rhoe0 / rho0;
        T0 = 1.e8;  // initial guess
        nse_T_abar_from_e(rho0, e0, Ye0, T0, abar);
    }

    // compute the plasma neutrino losses at t0

    amrex::Real snu{0.0};
    amrex::Real dsnudt{0.0};
    amrex::Real dsnudd{0.0};
    amrex::Real dsnuda{0.0};
    amrex::Real dsnudz{0.0};

    amrex::Real zbar = abar * Ye0;

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

    amrex::Real abar0_out = nse_state.abar;
    amrex::Real bea0_out = nse_state.bea;
    amrex::Real dyedt0 = nse_state.dyedt;

    // construct initial sources at t0

    amrex::Real rhoe_source{};
    if (integrator_rp::nse_include_enu_weak == 1) {
        rhoe_source = -rho0 * (nse_state.e_nu + snu);
    } else {
        rhoe_source = -rho0 * snu;
    }
    rhoe_source += rho0 * nse_state.dbeadt * C::MeV2eV * C::ev2erg * C::n_A;
    rhoe_source += C::n_A * (C::m_n - (C::m_p + C::m_e)) * C::c_light * C::c_light * rho0 * dyedt0;

    amrex::Real rhoaux_source[NumAux];
    rhoaux_source[iye] = rho0 * dyedt0;
    rhoaux_source[iabar] = 0.0;
    rhoaux_source[ibea] = rho0 * nse_state.dbeadt;

    // evolve for eps * dt;

    amrex::Real tau = integrator_rp::nse_deriv_dt_factor * dt;

    amrex::Real rho1 = rho0 + tau * ydot_a[SRHO];
    amrex::Real rhoe1 = rhoe0 + tau * ydot_a[SEINT] + tau * rhoe_source;
    amrex::Real rhoaux1[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        rhoaux1[n] = rhoaux0[n] + tau * ydot_a[SFX+n] + tau * rhoaux_source[n];
    }

    // compute the temperature at t0 + tau

    amrex::Real T1;
    amrex::Real Ye1 = rhoaux1[iye] / rho1;

    if (T_fixed > 0) {
        T1 = T_fixed;
    } else {
        amrex::Real abar1_out{};
        amrex::Real e1 = rhoe1 / rho1;
        T1 = T0;
        nse_T_abar_from_e(rho1, e1, Ye1, T1, abar1_out);
    }

    // call NSE at t0 + tau

    nse_state.T = T1;
    nse_state.rho = rho1;
    nse_state.Ye = Ye1;

    nse_interp(nse_state, skip_X_fill);

    amrex::Real bea1_out = nse_state.bea;

    // update abar -- this will be the same as nse_T_abar_from_e, but
    // for the case with T_fixed > 0, this abar will be consistent
    // with NSE.
    amrex::Real abar1_out = nse_state.abar;

    // construct the finite-difference approximation to the derivatives

    // note that abar, (B/A) here have already seen advection
    // implicitly

    amrex::Real rho_bea_tilde = rho1 * bea1_out - tau * ydot_a[SFX+ibea];
    amrex::Real rho_dBEA = rho_bea_tilde - rho0 * bea0_out; // this is MeV / nucleon * g / cm**3

    amrex::Real rho_abar_tilde = rho1 * abar1_out - tau * ydot_a[SFX+iabar];
    amrex::Real rho_dabar = rho_abar_tilde - rho0 * abar0_out; // this is MeV / nucleon * g / cm**3

    drhoedt = rho_dBEA * C::MeV2eV * C::ev2erg * C::n_A / tau;
    drhoedt += C::n_A * (C::m_n - (C::m_p + C::m_e)) * C::c_light * C::c_light * rho0 * dyedt0;
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
void sdc_nse_burn(BurnT& state, const amrex::Real dt) {

    using namespace AuxZero;

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // store the initial state

    amrex::Real rho_old = state.y[SRHO];
    amrex::Real rhoe_old = state.y[SEINT];
    amrex::Real rhoaux_old[NumAux];
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

    amrex::Real drhoedt;
    amrex::Real drhoauxdt[NumAux];

    nse_derivs(rho_old, rhoe_old, rhoaux_old,
               dt, state.ydot_a,
               drhoedt, drhoauxdt, state.T_fixed);

    // evolve to the midpoint in time

    amrex::Real rho_tmp = rho_old + 0.5_rt * dt * state.ydot_a[SRHO];
    amrex::Real rhoe_tmp = rhoe_old + 0.5_rt * dt * (state.ydot_a[SEINT] + drhoedt);
    amrex::Real rhoaux_tmp[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        rhoaux_tmp[n] = rhoaux_old[n] +
            0.5_rt * dt * (state.ydot_a[SFX+n] + drhoauxdt[n]);
    }

    // compute the derivatives at the midpoint in time

    nse_derivs(rho_tmp, rhoe_tmp, rhoaux_tmp,
               dt, state.ydot_a,
               drhoedt, drhoauxdt, state.T_fixed);

    // evolve to the new time

    amrex::Real rho_new = rho_old + dt * state.ydot_a[SRHO];
    amrex::Real rhoe_new = rhoe_old + dt * (state.ydot_a[SEINT] + drhoedt);
    amrex::Real rhoaux_new[NumAux];
    for (int n = 0; n < NumAux; ++n) {
        rhoaux_new[n] = rhoaux_old[n] +
            dt * (state.ydot_a[SFX+n] + drhoauxdt[n]);
    }

    // get the new temperature.  This will also get a new abar, but
    // we'll use the version from the NSE table compute below for the
    // final abar.

    amrex::Real T_new;
    amrex::Real abar_new;
    amrex::Real Ye_new = rhoaux_new[iye] / rho_new;

    if (state.T_fixed > 0) {
        T_new = state.T_fixed;
        // note: this is not consistent with NSE
        abar_new = rhoaux_new[iabar] / rho_new;
    } else {
        amrex::Real e_new = rhoe_new / rho_new;
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

    amrex::Real sum_X{0.0_rt};
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

    amrex::Real rho_enucdot = (rhoe_new - rhoe_old) / dt - state.ydot_a[SEINT];

    state.y[SEINT] = rhoe_new;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_enucdot;

}

#else
/// Self-Consistent NSE VERSION

///
/// this acts as an explicit Euler step for the system (rho e, rho)
/// on input, *_source are the reactive sources at time t0 and on output
/// they are the sources at time t0+dt
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_derivs(const amrex::Real rho0, const amrex::Real rhoe0,
                const amrex::Real rhoYe0, const amrex::Real T_old,
                const amrex::Real dt, amrex::Real &mu_p, amrex::Real &mu_n,
                amrex::Real &drhoyedt_weak, const amrex::Real drhoyedt_a,
                amrex::Array1D<amrex::Real, 1, neqs> &ydot_weak,
                const amrex::Real *ydot_a, amrex::Real& drhoedt,
                const amrex::Real T_fixed) {

    // start with the current state and compute T

    amrex::Real T0;
    amrex::Real Ye0 = rhoYe0 / rho0;
    amrex::Real abar{0.0_rt};

    // Get initial temperature and abar for

    if (T_fixed > 0) {
        T0 = T_fixed;
    } else {
        // initial guess for eos
        T0 = T_old;
        amrex::Real e0 = rhoe0 / rho0;
        nse_T_abar_from_e(rho0, e0, Ye0, T0, abar, mu_p, mu_n);
    }

    // initialize burn_state for NSE state

    burn_t burn_state;

    burn_state.T = T0;
    burn_state.rho = rho0;
    burn_state.y_e = Ye0;
    burn_state.mu_p = mu_p;
    burn_state.mu_n = mu_n;

    // get NSE state at t0

    auto nse_state0 = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

    // Still need to get abar for T_fixed > 0

    if (T_fixed > 0) {
        for (int n; n < NumSpec; ++n) {
            abar += nse_state0.xn[n] * aion_inv[n];
        }
        abar = 1.0_rt / abar;
    }

    // compute the plasma neutrino losses
    // here abar should already be in-sync with NSE

    amrex::Real snu{0.0_rt};
    amrex::Real dsnudt{0.0_rt};
    amrex::Real dsnudd{0.0_rt};
    amrex::Real dsnuda{0.0_rt};
    amrex::Real dsnudz{0.0_rt};

    amrex::Real zbar = abar * Ye0;

#ifdef NEUTRINOS
    {
        constexpr int do_derivatives = 0;
        sneut5<do_derivatives>(T0, rho0, abar, zbar,
                               snu, dsnudt, dsnudd,
                               dsnuda, dsnudz);
    }
#endif

    // Get molar fractions, ydots and neutrino loss due to weak rates

    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    amrex::Real e_nu {0.0_rt};

    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = nse_state0.xn[n-1] * aion_inv[n-1];
        // Y(n) = Xn0[n-1] * aion_inv[n-1];
    }

    // Fill in ydot with only weak rates contributing

    get_ydot_weak(nse_state0, ydot_weak, e_nu, Y);

    // rate_t rate_eval;
    // constexpr int do_T_derivatives = 0;
    // evaluate_rates<do_T_derivatives, rate_t>(burn_state, rate_eval);
    // rhs_nuc(burn_state, ydot_weak, Y, rate_eval.screened_rates);
    // e_nu = rate_eval.enuc_weak;


    ///
    /// construct initial sources at t0
    ///

    // Find drhoyedt_weak from weak reaction

    drhoyedt_weak = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        drhoyedt_weak += rho0 * zion[n] * ydot_weak(n+1);
    }

    // find rhoe_source

    amrex::Real rhoe_source{0.0_rt};
    for (int n = 1; n <= NumSpec; ++n) {
        rhoe_source += rho0 * ydot_weak(n) * network::mion(n);
    }
    rhoe_source *= C::Legacy::enuc_conv2;

    // include plasma neutrino terms

    rhoe_source -= rho0 * snu;

    // include neutrino loss terms from weak rates
    // Note that e_nu here is already negative so we add here.

    if (integrator_rp::nse_include_enu_weak == 1) {
        rhoe_source += rho0 * e_nu;
    }

    //
    // evolve for eps * dt;
    //

    amrex::Real tau = integrator_rp::nse_deriv_dt_factor * dt;

    amrex::Real rho1 = rho0 + tau * ydot_a[SRHO];
    amrex::Real rhoe1 = rhoe0 + tau * (ydot_a[SEINT] + rhoe_source);
    amrex::Real rhoYe1 = rhoYe0 + tau * (drhoyedt_weak + drhoyedt_a);

    // compute the temperature at t0 + tau

    amrex::Real T1;
    amrex::Real Ye1 = rhoYe1 / rho1;

    if (T_fixed > 0) {
        T1 = T_fixed;
    } else {
        // initial guess for eos
        T1 = T0;
        amrex::Real e1 = rhoe1 / rho1;
        amrex::Real abar1;
        nse_T_abar_from_e(rho1, e1, Ye1, T1, abar1, mu_p, mu_n);
    }

    // update burn_state at t0 + tau

    burn_state.T = T1;
    burn_state.rho = rho1;
    burn_state.y_e = Ye1;

    // call NSE at t0 + tau with

    // Note this new NSE state has advection contribution

    auto nse_state1 = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

    // construct the finite-difference approximation to the derivatives

    // subtract off advection contribution to get pure drhoedt from reaction

    drhoedt = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        drhoedt += (rho1 * nse_state1.xn[n] - tau * ydot_a[SFS+n] -
                    rho0 * nse_state0.xn[n]) * aion_inv[n] * network::mion(n+1) / tau;
        // drhoedt += (rho1 * nse_state1.xn[n] - tau * ydot_a[SFS+n] -
        //             rho0 * Xn0[n]) * aion_inv[n] * network::mion(n+1) / tau;
    }
    drhoedt *= C::Legacy::enuc_conv2;

    // include neutrino terms again

    drhoedt -= rho1 * snu;
    if (integrator_rp::nse_include_enu_weak == 1) {
        drhoedt += rho1 * e_nu;
    }

}


///
/// update the state due to NSE changes for the simplified-SDC case
/// this version works with the self-consistent NSE
///

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void sdc_nse_burn(BurnT& state, const amrex::Real dt) {

    state.success = true;
    state.n_rhs = 0;
    state.n_jac = 0;

    // store the initial state

    amrex::Real rho_old = state.y[SRHO];
    amrex::Real rhoe_old = state.y[SEINT];
    amrex::Real rhoye_old = rho_old * state.y_e;
    amrex::Real T_old = state.T;
    amrex::Real mu_p = state.mu_p;
    amrex::Real mu_n = state.mu_n;

    // density and momentum have no reactive sources

    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    // do an RK2 integration

    // get the derivatives, rho_dyedt, drhoedt, ydot_weak at t = t^n

    amrex::Array1D<amrex::Real, 1, neqs> ydot_weak;
    amrex::Real drhoedt{0.0_rt};
    amrex::Real drhoyedt_weak{0.0_rt};
    amrex::Real drhoyedt_a{0.0_rt};

    // find the advection contribution: drhoyedt_a and drhoabardt_a

    for (int n = 0; n < NumSpec; ++n) {
        drhoyedt_a += state.ydot_a[SFS+n] * zion[n] * aion_inv[n];
    }

    nse_derivs(rho_old, rhoe_old,
               rhoye_old, T_old,
               dt, mu_p, mu_n,
               drhoyedt_weak, drhoyedt_a,
               ydot_weak, state.ydot_a,
               drhoedt, state.T_fixed);

    // evolve to the midpoint in time

    amrex::Real rho_tmp = rho_old + 0.5_rt * dt * state.ydot_a[SRHO];
    amrex::Real rhoe_tmp = rhoe_old + 0.5_rt * dt * (state.ydot_a[SEINT] + drhoedt);
    amrex::Real rhoye_tmp = rhoye_old + 0.5_rt * dt * (drhoyedt_weak + drhoyedt_a);

    // compute the derivatives at the midpoint in time

    nse_derivs(rho_tmp, rhoe_tmp,
               rhoye_tmp, T_old,
               dt, mu_p, mu_n,
               drhoyedt_weak, drhoyedt_a,
               ydot_weak, state.ydot_a,
               drhoedt, state.T_fixed);

    // evolve to the new time

    amrex::Real rho_new = rho_old + dt * state.ydot_a[SRHO];
    amrex::Real rhoe_new = rhoe_old + dt * (state.ydot_a[SEINT] + drhoedt);
    amrex::Real rhoye_new = rhoye_old + dt * (drhoyedt_weak + drhoyedt_a);

    // Update the temperature with new internal energy

    amrex::Real T_new;
    amrex::Real Ye_new = rhoye_new / rho_new;

    if (state.T_fixed > 0) {
        T_new = state.T_fixed;
    } else {
        // initial eos guess
        T_new = T_old;
        amrex::Real e_new = rhoe_new / rho_new;
        amrex::Real abar_new;
        nse_T_abar_from_e(rho_new, e_new, Ye_new, T_new, abar_new, mu_p, mu_n);
    }

    // With updated temp, rho, and Ye get final NSE state

    burn_t burn_state;
    burn_state.T = T_new;
    burn_state.rho = rho_new;
    burn_state.y_e = Ye_new;
    burn_state.mu_p = mu_p;
    burn_state.mu_n = mu_n;

    auto nse_state = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

    // store the new state
    // make sure they are normalized

    amrex::Real sum_X{0.0_rt};
    for (auto & xn : nse_state.xn) {
        xn = std::clamp(xn, small_x, 1.0_rt);
        sum_X += xn;
    }

    for (auto & xn : nse_state.xn) {
        xn /= sum_X;
    }

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = rho_new * nse_state.xn[n];
    }

    // density and momenta have already been updated

    // get the energy release from the change in rhoe over the timestep
    // excluding any advection, and use that to update the total energy

    amrex::Real rho_enucdot = (rhoe_new - rhoe_old) / dt - state.ydot_a[SEINT];

    state.y[SEINT] = rhoe_new;
    state.y[SEDEN] += dt * state.ydot_a[SEDEN] + dt * rho_enucdot;

    // store the chemical potentials

    state.mu_p = nse_state.mu_p;
    state.mu_n = nse_state.mu_n;
}
#endif

#endif
