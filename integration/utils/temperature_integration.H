#ifndef _temperature_integration_H_
#define _temperature_integration_H_

#include <network.H>
#include <burn_type.H>
#include <extern_parameters.H>

// Sets up the temperature equation. This should be called from
// within the actual_rhs routine but is provided here as a convenience
// since most networks will use the same temperature ODE.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void temperature_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    // Set up the temperature ODE.  For constant pressure, Dp/Dt = 0, we
    // evolve :
    //    dT/dt = (1/c_p) [ Hnuc ]
    //
    // For constant volume, div{U} = 0, and we evolve:
    //    dT/dt = (1/c_v) [ Hnuc ]
    //
    // See low Mach paper III, including Eq. A3 for details.
    // Note that we no longer include the chemical potential (dE/dX or dH/dX)
    // terms because we believe they analytically should vanish.

    if (state.self_heat) {

        if (do_constant_volume_burn) {

            Real cv;

            if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

                cv = state.cv + (state.T - state.T_old) * state.dcvdT;

            }
            else {

                cv = state.cv;

            }

            Real cvInv = 1.0_rt / cv;

            ydot(net_itemp) = ydot(net_ienuc) * cvInv;

        }
        else {

            Real cp;

            if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

                cp = state.cp + (state.T - state.T_old) * state.dcpdT;

            }
            else {

                cp = state.cp;

            }

            Real cpInv = 1.0_rt / cp;

            ydot(net_itemp) = ydot(net_ienuc) * cpInv;

        }

    }

}


// Sets up the temperature entries in the Jacobian. This should be called from
// within the actual_jac routine but is provided here as a convenience
// since most networks will use the same temperature ODE.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void temperature_jac (burn_t& state, Array2D<Real, 1, neqs, 1, neqs>& jac)
{

    if (state.self_heat) {

        Real cspec, cspecInv;

        // Temperature Jacobian elements

        if (do_constant_volume_burn) {

            if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

                cspec = state.cv + (state.T - state.T_old) * state.dcvdT;

            }
            else {

                cspec = state.cv;

            }

        }
        else {

            if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

                cspec = state.cp + (state.T - state.T_old) * state.dcpdT;

            }
            else {

                cspec = state.cp;

            }

        }

        cspecInv = 1.0_rt / cspec;

        // d(itemp) / d(yi)

        for (int k = 1; k <= NumSpec; ++k) {
            jac(net_itemp, k) = jac(net_ienuc, k) * cspecInv;
        }

        // d(itemp) / d(temp) -- we get this from the equation for d (denuc / dt) / dT
        // since dT / dt = 1 / c_x denuc / dt in our formalism

        jac(net_itemp, net_itemp) = jac(net_ienuc, net_itemp) * cspecInv;

        // d(itemp) / d(enuc)

        jac(net_itemp, net_ienuc) = 0.0_rt;

    }

}

#endif
