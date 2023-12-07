#ifndef NSE_EOS_H
#define NSE_EOS_H

#include <AMReX_REAL.H>

#include <eos.H>

#include <extern_parameters.H>
#include <nse_table_type.H>
#include <nse_table.H>


using namespace amrex;


///
/// if we are in NSE, then the entire thermodynamic state is just
/// a function of rho, T, Ye.  We can write the energy as:
///
///    e = e(rho, T, Y_e, Abar(rho, T, Ye))
///
/// where we note that Abar is a function of those same inputs.
/// This function inverts this form of the EOS to find the T
/// that satisfies the EOS and NSE.
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real
nse_T_from_e(const Real rho, const Real e_in, const Real Ye, const Real T_guess) {

    using namespace AuxZero;

    const Real ttol{1.e-8_rt};
    const int max_iter{100};
    const Real eps{1.e-8_rt};

    // we need the full EOS type, since we need de/dA
    //eos_extra_t eos_state;

    bool converged{false};

    Real T{T_guess};
    int iter{};

    while (not converged && iter < max_iter) {

        // call NSE table to get Abar
        nse_table_t nse_state;
        nse_state.T = T;
        nse_state.rho = rho;
        nse_state.Ye = Ye;

        constexpr bool skip_X_fill{true};
        nse_interp(nse_state, skip_X_fill);
        Real abar_old = nse_state.abar;

        // call the EOS with the initial guess for T

        eos_state.rho = rho;
        eos_state.T = T;
        eos_state.aux[iye] = Ye;
        eos_state.aux[iabar] = abar_old;
        eos(eos_input_rt, eos_state);

        // f is the quantity we want to zero

        Real f = eos_state.e - e_in;

        // perturb T and call NSE again to get perturbed Abar
        nse_state.T *= (1.0_rt + eps);
        nse_interp(nse_state,skip_X_fill);

        // estimate dAbar/dT
        Real dabar_dT = (nse_state.abar - abar_old) / (eps * T);

        // compute the correction to our guess
        Real dT = f / (eos_state.dedT + eos_state.dedA * dabar_dT);

        // update the temperature
        T = std::clamp(T + dT, 0.5 * T, 2.0 * T);

        // check convergence

        if (std::abs(dT) < ttol * T) {
            converged = true;
        }
        iter++;
    }

    return T;

}

#endif
