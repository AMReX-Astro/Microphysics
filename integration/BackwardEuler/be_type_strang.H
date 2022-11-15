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
void evaluate_rhs(burn_t state, Array1D<Real, 1, BE_NEQS>& ydot)
{

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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real initial_dt (burn_t& state, const Real tend, Array1D<Real, 1, BE_NEQS>& ydot)
{
    // this is a version of the initial timestep estimator from VODE

    Array1D<Real, 1, BE_NEQS> y;

    for (int n = 1; n <= NumSpec; n++) {
        y(n) = state.xn[n-1];
    }
    y(net_ienuc) = state.e;

    // initial lower and upper bounds on the timestep

    Real hL = 100.0_rt * std::numeric_limits<Real>::epsilon() * tend;
    Real hU = 0.1_rt * tend;

    // initial guess for the itteraction

    Real h = std::sqrt(hL * hU);
    Real h_old = 10.0_rt * h;

    // Iterate on ddydtt = (RHS(t + h, y + h * dydt) - dydt) / h

    Array1D<Real, 1, BE_NEQS> ewt;
    Array1D<Real, 1, BE_NEQS> y_temp;
    Array1D<Real, 1, BE_NEQS> ydot_temp;
    Array1D<Real, 1, BE_NEQS> ddydtt;

    for (int n = 1; n <= 4; n++) {

       h_old = h;

       // Get the error weighting -- this is similar to VODE's dewset
       // routine

       for (int ii = 1; ii <= NumSpec; ii++) {
           ewt(ii) = rtol_spec * std::abs(y(ii)) + atol_spec;
       }
       ewt(net_ienuc) = rtol_enuc * std::abs(y(net_ienuc)) + atol_enuc;

       // Construct the trial point.

       for (int ii = 1; ii <= BE_NEQS; ii++) {
           y_temp(ii) = y(ii) + h * ydot(ii);
       }

       // Call the RHS, then estimate the finite difference.

       burn_t state_temp = state;
       for (int ii = 1; ii <= NumSpec; ii++) {
           state_temp.xn[ii-1] = y_temp(ii);
       }
       state_temp.e = y_temp(net_ienuc);

       clean_state(state_temp);

       evaluate_rhs(state_temp, ydot_temp);

       for (int ii = 1; ii <= BE_NEQS; ii++) {
           ddydtt(ii) = (ydot_temp(ii) - ydot(ii)) / h;
       }

       Real yddnorm = 0.0_rt;
       for (int ii = 1; ii <= BE_NEQS; ii++) {
           yddnorm += std::pow(ddydtt(ii) * ewt(ii), 2);
       }
       yddnorm = std::sqrt(yddnorm / BE_NEQS);

       if (yddnorm*hU*hU > 2.0_rt) {
          h = std::sqrt(2.0_rt / yddnorm);
       } else {
          h = std::sqrt(h * hU);
       }

       if (h_old < 2.0_rt * h && h_old > 0.5_rt * h) {
           break;
       }

    }

    // Save the final timestep, with a bias factor.

    Real dt = h / 2.0_rt;
    dt = amrex::min(amrex::max(h, hL), hU);

    dt = amrex::min(dt, ode_max_dt);

    return dt;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void initialize_state (burn_t& state)
{
    // We assume that (rho, T) coming in are valid, do an EOS call
    // to fill the rest of the thermodynamic variables.

    eos(eos_input_rt, state);

    state.success = true;

    state.n_rhs = 0;
}

#endif
