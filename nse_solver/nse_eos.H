#ifndef NSE_EOS_H
#define NSE_EOS_H

#include <AMReX_REAL.H>

#include <eos.H>

#include <extern_parameters.H>
#include <nse_solver.H>

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
nse_T_abar_from_e(const Real rho, const Real e_in, const Real Ye,
                  Real& T, Real& abar) {

    constexpr Real ttol{1.e-6_rt};
    constexpr int max_iter{100};

    // we need the full EOS type, since we need de/dA
    //eos_extra_t eos_state;

    bool converged{false};

    int iter{0};

    while (not converged && iter < max_iter) {

        // fill in Temp, rho, and Ye to get NSE state

        burn_t burn_state;

        burn_state.T = T;
        burn_state.rho = rho;
        burn_state.Ye = Ye;

        auto nse_state = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

        // amrex::Real abar_old{0.0_rt];

        // for (int n = 0; n < NumSpec; ++n) {
        //     abar_old = nse_state.xn[n] * aion_inv[n];
        // }
        // abar_old = 1.0_rt / abar_old;

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

        amrex::Real dabar_dT = nse_dabar_dT(T, rho, Ye);

        // compute the correction to our guess

        Real dT = -f / (eos_state.dedT + eos_state.dedA * dabar_dT
                                  + Ye * eos_state.dedZ * dabar_dT);

        // update the temperature and abar

        T = std::clamp(T + dT, 0.25 * T, 4.0 * T);
        abar = eos_state.abar;

        // check convergence

        if (std::abs(dT) < ttol * T) {
            converged = true;
        }
        iter++;

    }

}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real nse_abar(const amrex::Real T,
                     const amrex::Real rho,
                     const amrex::Real Ye) {
    ///
    /// This function calculates abar from NSE using
    /// Temp, rho, and Ye
    ///

    burn_t burn_state;
    burn_state.rho = rho;
    burn_state.y_e = Ye;
    burn_state.T = T;

    auto nse_state = get_actual_nse_state(burn_state, 1.0e-10_rt, true);

    amrex::Real abar{0.0_rt}
    for (int n = 0; n < NumSpec; ++n) {
        abar += nse_state.xn[n] * aion_inv[n];
    }

    return abar;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::Real nse_dabar_dT(const amrex::Real T,
                         const amrex::Real rho,
                         const amrex::Real Ye) {
    ///
    /// This function calculates dabar_dT using
    /// five-point stencil method.
    /// This should be 4th order accurate.
    ///

    // deviation in temperature

    const amrex::Real dT = 0.001_rt * T;

    // Caclulate derivative using five-point stencil method

    amrex::Real dabar_dT = (-nse_abar(T + 2.0_rt*dT, rho, Ye) +
                             8.0_rt * nse_abar(T+dT, rho, Ye) -
                             8.0_rt * nse_abar(T-dT, rho, Ye) +
                             nse_abar(T - 2.0_rt*dT, rho, Ye)
                            ) / (12.0_rt * dT);

    return dabar_dT
}

#endif
