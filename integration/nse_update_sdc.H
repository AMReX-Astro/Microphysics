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
#include <neutrino.H>

#ifdef NSE_TABLE
#include <nse_table.H>
#include <nse_eos.H>
#endif
#ifdef NSE_NET
#include <nse_solver.H>
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
    amrex::Real pair{0.0};
    amrex::Real phot{0.0};
    amrex::Real plas{0.0};
    amrex::Real brem{0.0};

    amrex::Real zbar = abar * Ye0;

#ifdef NEUTRINOS
    constexpr int do_derivatives = 0;
    neutrino_cooling<do_derivatives>(T0, rho0, abar, zbar,
                                     snu, dsnudt, dsnudd, dsnuda, dsnudz,
                                     pair, phot, plas, brem);
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
        xn = amrex::Clamp(xn, small_x, 1.0_rt);
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
/// This computes drhoedt and drhoyedt_weak for a given (rho, rho e, and rho Ye)
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_derivs(const amrex::Real rho0, const amrex::Real rhoe0,
                const amrex::Real rhoYe0, amrex::Real &T0,
                amrex::Real &mu_p, amrex::Real &mu_n,
                amrex::Real &drhoyedt_weak, amrex::Real& drhoedt,
                const amrex::Real T_fixed) {

    // initialize burn_state for NSE state

    amrex::Real Ye0 = rhoYe0 / rho0;

    burn_t burn_state;

    burn_state.T = (T_fixed > 0.0_rt) ? T_fixed : T0;
    burn_state.rho = rho0;
    burn_state.e = rhoe0 / rho0;
    burn_state.y_e = Ye0;
    burn_state.mu_p = mu_p;
    burn_state.mu_n = mu_n;

    // Since internal energy is updated during integration,
    // we use (rho, e, Ye) as input to find the NSE state.
    // If T is fixed, i.e. T_fixed > 0, then use (rho, T, Ye)

    auto nse_input = (T_fixed > 0.0_rt) ? nse_input_rty : nse_input_rey;
    auto nse_state = get_actual_nse_state(nse_input, burn_state,
                                          1.0e-10_rt, true);

    // Update the input temperature, mu_p and mu_n with the new NSE state.

    if (T_fixed <= 0.0_rt) {
        T0 = nse_state.T;
    }
    mu_p = nse_state.mu_p;
    mu_n = nse_state.mu_n;

#ifdef NEUTRINOS
    // compute the plasma neutrino losses

    // Get abar for neutrino

    amrex::Real abar{0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        abar += nse_state.xn[n] * aion_inv[n];
    }
    abar = 1.0_rt / abar;

    amrex::Real snu{0.0_rt};
    amrex::Real dsnudt{0.0_rt};
    amrex::Real dsnudd{0.0_rt};
    amrex::Real dsnuda{0.0_rt};
    amrex::Real dsnudz{0.0_rt};
    amrex::Real pair{0.0};
    amrex::Real phot{0.0};
    amrex::Real plas{0.0};
    amrex::Real brem{0.0};

    amrex::Real zbar = abar * Ye0;

    constexpr int do_derivatives = 0;
    neutrino_cooling<do_derivatives>(T0, rho0, abar, zbar,
                                     snu, dsnudt, dsnudd, dsnuda, dsnudz,
                                     pair, phot, plas, brem);
#endif

    // Get molar fractions, ydots and neutrino loss due to weak rates

    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    amrex::Real e_nu {0.0_rt};

    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = nse_state.xn[n-1] * aion_inv[n-1];
    }

    // Fill in ydot with only weak rates contributing

    amrex::Array1D<amrex::Real, 1, neqs> ydot_weak;
    get_ydot_weak(nse_state, ydot_weak, e_nu, Y);

    // Find drhoyedt_weak from weak reaction

    drhoyedt_weak = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        drhoyedt_weak += rho0 * zion[n] * ydot_weak(n+1);
    }

    // Find drhoedt

    drhoedt = 0.0_rt;
    for (int n = 1; n <= NumSpec; ++n) {
        drhoedt += rho0 * ydot_weak(n) * network::mion(n);
    }
    drhoedt *= C::enuc_conv2;

    // include plasma neutrino terms

    drhoedt -= rho0 * snu;

    // include neutrino loss terms from weak rates
    // Note that e_nu here is already negative so we add here.

    if (integrator_rp::nse_include_enu_weak == 1) {
        drhoedt += rho0 * e_nu;
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

    // do an RK2 integration

    // initialize derivatives needed

    amrex::Real drhoedt{0.0_rt};
    amrex::Real drhoyedt_weak{0.0_rt};
    amrex::Real drhoyedt_a{0.0_rt};

    // find the advection contribution for ye: drhoyedt_a
    // The advection contribution is already time-centered

    for (int n = 0; n < NumSpec; ++n) {
        drhoyedt_a += state.ydot_a[SFS+n] * zion[n] * aion_inv[n];
    }

    // get the derivatives: drhoyedt_weak and drhoedt at t = t^n

    nse_derivs(rho_old, rhoe_old,
               rhoye_old, T_old,
               mu_p, mu_n,
               drhoyedt_weak, drhoedt,
               state.T_fixed);

    // evolve to the midpoint in time and the midpoint NSE state.

    amrex::Real rho_tmp = rho_old + 0.5_rt * dt * state.ydot_a[SRHO];
    amrex::Real rhoe_tmp = rhoe_old + 0.5_rt * dt * (state.ydot_a[SEINT] + drhoedt);
    amrex::Real rhoye_tmp = rhoye_old + 0.5_rt * dt * (drhoyedt_weak + drhoyedt_a);

    // compute the derivatives at the midpoint in time

    nse_derivs(rho_tmp, rhoe_tmp,
               rhoye_tmp, T_old,
               mu_p, mu_n,
               drhoyedt_weak, drhoedt,
               state.T_fixed);

    // evolve to the new time

    amrex::Real rho_new = rho_old + dt * state.ydot_a[SRHO];
    amrex::Real rhoe_new = rhoe_old + dt * (state.ydot_a[SEINT] + drhoedt);
    amrex::Real rhoye_new = rhoye_old + dt * (drhoyedt_weak + drhoyedt_a);

    // With updated temp, rho, and Ye get final NSE state

    burn_t burn_state;
    burn_state.T = (state.T_fixed > 0.0_rt) ? state.T_fixed : T_old;
    burn_state.rho = rho_new;
    burn_state.e = rhoe_new / rho_new;
    burn_state.y_e = rhoye_new / rho_new;
    burn_state.mu_p = mu_p;
    burn_state.mu_n = mu_n;

    auto nse_input = (state.T_fixed > 0.0_rt) ? nse_input_rty : nse_input_rey;
    auto nse_state = get_actual_nse_state(nse_input, burn_state,
                                          1.0e-10_rt, true);

    // Update rho X using new NSE state
    // These should be already normalized due to NSE constraint

    for (int n = 0; n < NumSpec; ++n) {
        state.y[SFS+n] = nse_state.y[SFS+n];
    }

    // Update density and momentum, which have no reactive sources

    state.y[SRHO] += dt * state.ydot_a[SRHO];
    state.y[SMX] += dt * state.ydot_a[SMX];
    state.y[SMY] += dt * state.ydot_a[SMY];
    state.y[SMZ] += dt * state.ydot_a[SMZ];

    // Update total and internal energy.

    state.y[SEINT] = rhoe_new;
    state.y[SEDEN] += dt * (state.ydot_a[SEDEN] + drhoedt);

    // Update chemical potentials and temperature

    if (state.T_fixed <= 0.0_rt) {
        state.T = nse_state.T;
    }
    state.mu_p = nse_state.mu_p;
    state.mu_n = nse_state.mu_n;
}
#endif

#endif