#ifndef NSE_EOS_H
#define NSE_EOS_H

#include <AMReX_REAL.H>

#include <eos.H>

#include <extern_parameters.H>
#include <nse_table_type.H>
#include <nse_table.H>


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
/// The basic idea is that Abar and Zbar are both functions of
/// rho, T, Ye through the NSE table, so we express the energy
/// as:
///
///      e = e(rho, T, Abar(rho, T, Ye), Zbar(rho, T, Ye)
///
/// and NR on that.  Note that Zbar = Ye Abar, so we can group
/// those derivative terms together.
///
/// T and abar come in as initial guesses and are updated
/// on output
///
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
nse_T_abar_from_e(const Real rho, const Real e_in, const Real Ye,
                  Real& T, Real& abar) {

    using namespace amrex;
    using namespace AuxZero;

    const Real ttol{1.e-6_rt};
    const int max_iter{100};

    // we need the full EOS type, since we need de/dA
    //eos_extra_t eos_state;

    bool converged{false};

    int iter{};

    nse_table_t nse_state;

    while (not converged && iter < max_iter) {

        // call NSE table to get Abar
        nse_state.T = T;
        nse_state.rho = rho;
        nse_state.Ye = Ye;

        constexpr bool skip_X_fill{true};
        nse_interp(nse_state, skip_X_fill);
        Real abar_old = nse_state.abar;

        // call the EOS with the initial guess for T
        eos_extra_t eos_state;
        eos_state.rho = rho;
        eos_state.T = T;
        eos_state.aux[iye] = Ye;
        eos_state.aux[iabar] = abar_old;
        eos(eos_input_rt, eos_state);

        // f is the quantity we want to zero

        Real f = eos_state.e - e_in;

        Real dabar_dT = nse_interp_dT(T, rho, Ye, nse_table::abartab);

        // compute the correction to our guess

        Real dT = -f / (eos_state.dedT + eos_state.dedA * dabar_dT
                                  + Ye * eos_state.dedZ * dabar_dT);

        // update the temperature
        T = std::clamp(T + dT, 0.25 * T, 4.0 * T);

        // check convergence

        if (std::abs(dT) < ttol * T) {
            converged = true;
        }
        iter++;
    }

    // T is set to the last T
    // we just need to save abar for output
    abar = nse_state.abar;

}

#endif
