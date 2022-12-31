#ifndef actual_integrator_H
#define actual_integrator_H

#include <network.H>
#include <actual_network.H>
#include <rhs.H>
#include <actual_rhs.H>
#include <burn_type.H>
#include <rate_type.H>
#include <eos_type.H>
#include <eos.H>
#include <extern_parameters.H>

Real calculate_dt(burn_t& state, Array1D<Real, 1, NumSpec>& spec_rhs, Real& ener_rhs)
{
    //The idea of calculate_dt is to prevent a sharp change of dX/dt by resizing the time step.
    Real dt = 1.0e200_rt;

    for(int n=1; n<=NumSpec; ++n){

        if(state.xn[n-1] >= atol_spec){

            Real target_dX;

            if(spec_rhs(n) > 0.0_rt){
                // When dX_dt > 0, we limit dX, such that X + (dX/dt) * dt <= 1.0
                // The philosophy is the following: In general for creation terms (dX/dt > 0.0)
                // we have dX/dt ~ aX. .ie dX/dt is linear in X. Therefore we can limit dX
                target_dX = (maximum_timestep_change_factor - 1.0_rt) * state.xn[n-1];
            } else {
                target_dX = (1.0_rt - 1.0_rt / maximum_timestep_change_factor) * state.xn[n-1];
            }

            Real dX_dt = amrex::max(std::abs(spec_rhs(n)), 1.0e-30_rt);

            dt = amrex::min(dt, target_dX / dX_dt);
        }

    }

    if(integrate_energy){

        Real target_de;

        if(ener_rhs > 0.0){
            target_de = (maximum_timestep_change_factor - 1.0_rt) * state.e;
        } else {
            target_de = (1 - 1.0_rt / maximum_timestep_change_factor) * state.e;
        }

        Real de_dt = amrex::max(std::abs(ener_rhs), 1.0e-30_rt);

        dt = amrex::min(dt, target_de / de_dt);
    }
}


void initialize_state(burn_t& state){

    eos(eos_input_rt, state);

    state.sucess = true;

    state.n_rhs = 0;

}

void clean_state(burn_t& state){

    // We start by normalizing the abundances
    normalize_abundances_burn(state);

    //From the energy and density we compute the temperature,
    eos(eos_input_re, state);

    //And finally, we setup a limit to the temperature
    state.T = amrex::min(MAX_TEMP, amrex::max(state.T, EOSData::mintemp));
}

void evaluate_rhs(burn_t& state, Array1D<Real, 1, NumSpec>& spec_rhs, Real& ener_rhs){

    Array1D<Real, 1, neqs> ydot;

    actual_rhs(state, ydot);

    for(int n=1; n<=NumSpec; ++n){
        spec_rhs(n) = aion[n] * ydot(n);
    }

    ener_rhs = ydot(net_ienuc);
}

void actual_integrator(burn_t& state, Real dt)
{
    initialize_state(state);

    Real T_in = state.T;
    Real e_in = state.e;
    Real x_in[NumSpec];

    for(int n=1; n<= NumSpec; ++n){
        x_in = state.xn[n-1];
    }

    Real timestep_safety_factor = 1.0e-12;

    int num_timesteps = 0;

    while(t < (1.0_rt - timestep_safety_factor)*dt && num_timesteps < ode_max_steps){



    }

}
