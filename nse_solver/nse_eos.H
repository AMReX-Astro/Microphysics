#ifndef NSE_EOS_H
#define NSE_EOS_H

#include <AMReX_REAL.H>

#include <eos.H>

#include <extern_parameters.H>
#include <nse_solver.H>


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real nse_abar(const amrex::Real T, const amrex::Real rho,
                     const amrex::Real Ye, amrex::Real &mu_p,
                     amrex::Real &mu_n) {
    ///
    /// This function calculates abar from NSE using
    /// Temp, rho, and Ye
    ///

    burn_t burn_state;
    burn_state.rho = rho;
    burn_state.y_e = Ye;
    burn_state.T = T;
    burn_state.mu_p = mu_p;
    burn_state.mu_n = mu_n;

    auto nse_state = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

    amrex::Real abar{0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        abar += nse_state.xn[n] * aion_inv[n];
    }

    // update mu_p and mu_n

    mu_p = burn_state.mu_p;
    mu_n = burn_state.mu_n;

    return abar;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real nse_dabar_dT(const amrex::Real T, const amrex::Real rho,
                         const amrex::Real Ye, amrex::Real &mu_p,
                         amrex::Real &mu_n) {
    ///
    /// This function calculates dabar_dT using
    /// five-point stencil method.
    /// This should be 4th order accurate.
    ///

    // deviation in temperature

    const amrex::Real dT = 1.0e-5_rt* T;

    // Caclulate derivative using five-point stencil method

    amrex::Real dabar_dT = (-nse_abar(T + 2.0_rt*dT, rho, Ye, mu_p, mu_n) +
                             8.0_rt * nse_abar(T+dT, rho, Ye, mu_p, mu_n) -
                             8.0_rt * nse_abar(T-dT, rho, Ye, mu_p, mu_n) +
                             nse_abar(T - 2.0_rt*dT, rho, Ye, mu_p, mu_n)
                            ) / (12.0_rt * dT);

    return dabar_dT;
}

///
/// This function inverts this form of the EOS to find the T
/// that satisfies the EOS and NSE given an input e and rho.
///
/// if we are in NSE, then the entire thermodynamic state is just
/// a function of rho, T, Ye.  We can write the energy as:
///
///    e = e(rho, T, Y_e, Abar(rho, T, Ye))
///
/// where we note that Abar is a function of those same inputs.
///
/// The basic idea is that Abar and Zbar are both functions of
/// rho, T, Ye through NSE calculations, so we express the energy
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
nse_T_abar_from_e(const amrex::Real rho, const amrex::Real e_in,
                  const amrex::Real Ye,
                  amrex::Real &T, amrex::Real &abar,
                  amrex::Real &mu_p, amrex::Real &mu_n) {

    constexpr Real ttol{1.e-6_rt};
    constexpr int max_iter{100};

    // we need the full EOS type, since we need de/dA
    //eos_extra_t eos_state;

    bool converged{false};

    int iter{0};

    // initialize burn_state
    burn_t burn_state;
    burn_state.rho = rho;
    burn_state.y_e = Ye;
    burn_state.mu_p = mu_p;
    burn_state.mu_n = mu_n;

    while (not converged && iter < max_iter) {

        // update Temperature

        burn_state.T = T;

        auto nse_state = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

        // call the EOS with the initial guess for T
        // as well as density and NSE mass fractions

        eos_extra_t eos_state;
        for  (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = nse_state.xn[n];
        }

        // Call the EOS to get internal energy
        // abar, zbar, etc are calculated automatically in EOS

        eos_state.rho = rho;
        eos_state.T = T;
        eos(eos_input_rt, eos_state);

        // f is the quantity we want to zero

        amrex::Real f = eos_state.e - e_in;

        amrex::Real dabar_dT = nse_dabar_dT(T, rho, Ye,
                                            nse_state.mu_p, nse_state.mu_n);

        // compute the correction to our guess

        Real dT = -f / (eos_state.dedT + eos_state.dedA * dabar_dT
                                  + Ye * eos_state.dedZ * dabar_dT);

        // update the temperature and abar

        T = std::clamp(T + dT, 0.25 * T, 4.0 * T);
        abar = eos_state.abar;

        // update mu_p and mu_n

        mu_p = nse_state.mu_p;
        mu_n = nse_state.mu_n;

        // check convergence

        if (std::abs(dT) < ttol * T) {
            converged = true;
        }
        iter++;

    }

}

#endif
