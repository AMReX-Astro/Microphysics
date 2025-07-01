#ifndef NSE_SOLVER_H
#define NSE_SOLVER_H

#include <fundamental_constants.H>
#include <AMReX_REAL.H>
#include <eos_type.H>
#include <eos.H>
#include <network.H>
#include <burn_type.H>
#include <extern_parameters.H>
#include <cmath>
#include <AMReX_Array.H>
#include <actual_network.H>
#ifndef NEW_NETWORK_IMPLEMENTATION
#include <actual_rhs.H>
#endif
#include <eos_composition.H>
#include <microphysics_sort.H>
#include <hybrj.H>
#include <screen.H>
#include <cctype>
#include <algorithm>
#include <nse_eqns.H>

using namespace nse_rp;

// Two different input mode for NSE solver
// nse_input_rty: (rho, T, Ye)
// nse_input_rey: (rho, e, Ye)

enum nse_input_t : std::uint8_t {
    nse_input_rty = 0,
    nse_input_rey
};


template <typename T>
struct nse_solver_data
{
    T state;
    amrex::Array1D<amrex::Real, 1, NumSpec> u_c;
};


template <int neqs, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_hybrid_solver(nse_solver_data<T>& state_data,
                       amrex::Real eps=1.0e-10_rt) {
    // state is the nse_state from get_nonexponent_nse_state

    hybrj_t<neqs> hj;

    // we'll take x[1] = mu_p, x[2] = mu_n

    hj.xtol = eps;
    hj.mode = 1;

    // random flag number for evaluation in for loop;

    int flag = 0;

    // Fine-tune variables

    amrex::Real dx;
    bool is_pos_new;
    bool is_pos_old = false;

    amrex::Array1D<amrex::Real, 1, neqs> f;
    amrex::Array1D<amrex::Real, 1, neqs> outer_x;
    amrex::Array1D<amrex::Real, 1, neqs> inner_x;

    outer_x(1) = state_data.state.mu_p;
    outer_x(2) = state_data.state.mu_n;

    // for (int j = 1; j <= 2; ++j) {
    //     hj.diag(j) = 1.0_rt;
    // }

    // fine tuning initial guesses

    for (int i = 0; i < 20; ++i) {

        dx = 0.5_rt;
        inner_x(1) = outer_x(1);
        inner_x(2) = outer_x(2);

        for (int j = 0; j < 20; ++j) {

            hj.x(1) = inner_x(1);
            hj.x(2) = inner_x(2);

            hybrj<neqs, nse_solver_data<T>>(hj, state_data,
                                            nse_fcn<nse_solver_data<T>>,
                                            nse_jcn<nse_solver_data<T>>);

            nse_fcn(hj.x, f, state_data, flag);

            if (std::abs(f(1)) < eps && std::abs(f(2)) < eps) {

                state_data.state.mu_p = hj.x(1);
                state_data.state.mu_n = hj.x(2);
                return;
            }

            is_pos_new = f(1) > 0.0_rt && f(2) > 0.0_rt;

            if (is_pos_old != is_pos_new) {
                dx *= 0.8_rt;
            }

            if (is_pos_new) {
                inner_x(1) -= dx;
                inner_x(2) -= dx;
            }
            else{
                inner_x(1) += dx;
                inner_x(2) += dx;
            }

            is_pos_old = is_pos_new;
        }

        outer_x(1) -= 0.5_rt;

    }

    // if (hj.info != 1) {
    //     amrex::Error("failed to solve");
    // }
#ifndef AMREX_USE_GPU
    std::cout << "NSE solver failed with these conditions: " << std::endl;
    std::cout << "Temperature: " << state_data.state.T << std::endl;
    std::cout << "Density: " << state_data.state.rho << std::endl;
    std::cout << "Ye: " << state_data.state.y_e << std::endl;
    std::cout << "Initial mu_p and mu_n: " << state_data.state.mu_p
              << ", " << state_data.state.mu_n << std::endl;
#endif

    amrex::Error("failed to solve");
}

// A newton-raphson solver for finding nse state used for calibrating
// chemical potential of proton and neutron

template <int neqs, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_nr_solver(nse_solver_data<T>& state_data,
                   amrex::Real eps=1.0e-10_rt) {
    // state_data is the state_data after from
    // get_nonexponent_nse_state and compute_coulomb_contribution

    // whether nse solver converged or not

    bool converged = false;

    amrex::Array1D<amrex::Real, 1, neqs> f;
    amrex::Array2D<amrex::Real, 1, neqs, 1, neqs> jac;
    amrex::Array1D<amrex::Real, 1, neqs> x;
    int flag = 0;

    x(1) = state_data.state.mu_p;
    x(2) = state_data.state.mu_n;

    nse_jcn(x, jac, state_data, flag);
    nse_fcn(x, f, state_data, flag);

    // store determinant for finding inverse jac
    amrex::Real det;

    // store inverse jacobian
    amrex::Array2D<amrex::Real, 1, neqs, 1, neqs> inverse_jac;

    // difference in chemical potential of proton and neutron
    amrex::Real d_mu_p = std::numeric_limits<amrex::Real>::max();
    amrex::Real d_mu_n = std::numeric_limits<amrex::Real>::max();

    // begin newton-raphson
    for (int i = 0; i < max_nse_iters; ++i) {

        // check if current state fulfills constraint equation

        if (std::abs(d_mu_p) < eps * std::abs(x(1)) &&
            std::abs(d_mu_n) < eps * std::abs(x(2))) {
            converged = true;
            state_data.state.mu_p = x(1);
            state_data.state.mu_n = x(2);
            break;
        }

        // Find the max of the jacobian used for scaling determinant to prevent digit overflow

        auto scale_fac = amrex::max(jac(2,2), amrex::max(jac(2,1), amrex::max(jac(1,1), jac(1,2))));

        // if jacobians are small, then no need for scaling

        if (scale_fac < 1.0e150) {
            scale_fac = 1.0_rt;
        }

        // Specific inverse 2x2 matrix, perhaps can write a function for solving n systems of equations.

        det = jac(1, 1) / scale_fac * jac(2, 2) - jac(1, 2) / scale_fac * jac(2, 1);

        // check if determinant is 0

        if (det == 0.0_rt) {
            amrex::Error("Jacobian is a singular matrix! Try a different initial guess!");
        }

        // find inverse jacobian

        inverse_jac(1, 1) =  jac(2, 2) / scale_fac / det;
        inverse_jac(1, 2) = -jac(1, 2) / scale_fac / det;
        inverse_jac(2, 1) = -jac(2, 1) / scale_fac / det;
        inverse_jac(2, 2) =  jac(1, 1) / scale_fac / det;

        // find the difference

        d_mu_p = -(f(1) * inverse_jac(1, 1) + f(2) * inverse_jac(1, 2));
        d_mu_n = -(f(1) * inverse_jac(2, 1) + f(2) * inverse_jac(2, 2));

        // if diff goes beyond 1.0e3_rt, likely that its not making good progress..

        if (std::abs(d_mu_p) > 1.0e3_rt || std::abs(d_mu_n) > 1.0e3_rt) {
            amrex::Error("Not making good progress, breaking");
        }

        // update new solution

        x(1) += d_mu_p;
        x(2) += d_mu_n;

        // check whether solution results in nan

        if (std::isnan(x(1)) || std::isnan(x(2))) {
            amrex::Error("Nan encountered, likely due to overflow in digits or not making good progress");
        }

        // update constraint

        nse_jcn(x, jac, state_data, flag);
        nse_fcn(x, f, state_data, flag);
    }

    if (!converged) {
        amrex::Error("NSE solver failed to converge!");
    }
}


template <int neqs, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_solve(nse_solver_data<T>& state_data,
               const amrex::Real eps=1.0e-10_rt) {
    // Interface for calling the Hybrid-Powell or Netwon-Raphson Solver

    if (use_hybrid_solver) {
        nse_hybrid_solver<neqs, T>(state_data, eps);
    }
    else {
        bool singular_network = true;
        for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
            if (n == NSE_INDEX::H1_index) {
                continue;
            }
#endif
            if (zion[n] != aion[n] - zion[n]) {
                singular_network = false;
            }
        }

        if (singular_network) {
            amrex::Error("This network always results in singular jacobian matrix, thus can't find nse mass fraction using nr!");
        }

        nse_nr_solver<neqs, T>(state_data, eps);
    }
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_nse_state_rty(T& state,
                    const amrex::Real eps=1.0e-10_rt) {
    // Find the NSE state using (rho, T, Ye) as input.
    // This returns state_data.state with the NSE mass fraction and
    // solved chemical potentials.
    // It also updates chemical potentials to the input state.

    // Initializes the initial state

    nse_solver_data<T> state_data = {state, {0.0_rt}};

    // Since Temperature is fixed in this mode,
    // Get nse_state without the exponent term

    state_data.state = get_nonexponent_nse_state(state);

    // if use chabrier1998 screening, calculate the coulomb correction term

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    compute_coulomb_contribution(state_data.u_c, state);
#endif

    // invoke newton-raphson or hybrj to solve chemical potential of proton and neutron
    // which are the exponent part of the nse mass fraction calculation

    // For (rho, T, Ye) input, there are two constraint eqns.

    constexpr int nConstraintEqns = 2;

    nse_solve<nConstraintEqns, T>(state_data, eps);

    // Apply exponent for calculating nse mass fractions

    apply_nse_exponent(state_data.state, state_data.u_c);

    // update mu_n and mu_p to input state

    state.mu_p = state_data.state.mu_p;
    state.mu_n = state_data.state.mu_n;

    // update rhoX in the output nse_state

    for (int n = 0; n < NumSpec; ++n) {
        state_data.state.y[SFS+n] = state.rho * state_data.state.xn[n];
    }

    return state_data.state;
}


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

    auto nse_state = get_nse_state_rty(burn_state, 1.0e-10_rt);

    amrex::Real abar{0.0_rt};
    for (int n = 0; n < NumSpec; ++n) {
        abar += nse_state.xn[n] * aion_inv[n];
    }

    abar = 1.0_rt / abar;

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
    /// This function constructs dabar_dT
    /// This should be 2nd order accurate
    ///

    // deviation in temperature

    const amrex::Real dT = 1.0e-6_rt * T;

    // Calculate derivative using five-point stencil method

    // amrex::Real dabar_dT = (-nse_abar(T + 2.0_rt*dT, rho, Ye, mu_p, mu_n) +
    //                          8.0_rt * nse_abar(T+dT, rho, Ye, mu_p, mu_n) -
    //                          8.0_rt * nse_abar(T-dT, rho, Ye, mu_p, mu_n) +
    //                          nse_abar(T - 2.0_rt*dT, rho, Ye, mu_p, mu_n)
    //                         ) / (12.0_rt * dT);

    // Calculate derivative using central differencing

    amrex::Real dabar_dT = 0.5_rt * (nse_abar(T + dT, rho, Ye, mu_p, mu_n) -
                                     nse_abar(T - dT, rho, Ye, mu_p, mu_n)) / dT;

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
/// T come in as initial guesses and are updated on output
///

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
nse_T_from_e(const amrex::Real rho, const amrex::Real e_in,
             const amrex::Real Ye,
             amrex::Real &T,
             amrex::Real &mu_p, amrex::Real &mu_n) {

    constexpr amrex::Real ttol{1.e-8_rt};
    constexpr int max_iter{100};

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

        auto nse_state = get_nse_state_rty(burn_state, 1.0e-10_rt);

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

        amrex::Real dT = -f / (eos_state.dedT + eos_state.dedA * dabar_dT
                               + Ye * eos_state.dedZ * dabar_dT);

        // update the temperature

        T = amrex::Clamp(T + dT, 0.25 * T, 4.0 * T);

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


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_nse_state_rey(T& state,
                    const amrex::Real eps=1.0e-10_rt) {
    // Find the NSE state using (rho, e, Ye) as input.

    // When using (rho, e, Ye) as input
    // We can solve in two ways:
    // 1. Find the NSE state using (rho, T, Ye), then do NR iterations
    // on that to find the correct T that matches e.
    // 2. Find the NSE state directly using (rho, e, Ye),
    // where now we have an additional constraint equation on e.
    // So there are three constraint equations and a 3x3 Jacobian.

    // Get T0 such that it is consistent with e.
    // This also updates mu_p and mu_n of the input state.

    amrex::Real T0 = state.T;
    nse_T_from_e(state.rho, state.e, state.y_e,
                 T0, state.mu_p, state.mu_n);

    // Create a new state with the updated Temperature
    // Then find NSE state via (rho, T, Ye)

    burn_t burn_state;
    burn_state.T = T0;
    burn_state.rho = state.rho;
    burn_state.y_e = state.y_e;
    burn_state.e = state.e;
    burn_state.mu_p = state.mu_p;
    burn_state.mu_n = state.mu_n;

    auto nse_state = get_nse_state_rty(burn_state, eps);

    // update mu_n and mu_p to input state

    state.mu_p = nse_state.mu_p;
    state.mu_n = nse_state.mu_n;

    return nse_state;
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_actual_nse_state(nse_input_t input,
                       T& state, const amrex::Real eps=1.0e-10_rt,
                       bool input_ye_is_valid=false) {
    ///
    /// Main interface for getting the NSE state.
    /// There are two input modes:
    ///
    /// 1. (rho, T, Ye): nse_input_rty
    /// In this mode, state.rho, state.T and state.y_e should be defined.
    ///
    /// 2. (rho, e, Ye): nse_input_rey
    /// In this mode, state.rho, state.e and state.y_e should be defined.
    /// It uses state.T as the initial guess, so that should be defined too.
    ///
    /// In both mode, if input_ye_is_valid=false or state.y_e is out of range,
    /// then y_e is re-calculated based on state.y[SFS], i.e. rho X's.
    ///
    /// The input state will be updated with new chemical potentials
    /// The output state will have updated X's, rhoX's and
    /// Temperature, if using (rho, e, Ye) mode.
    ///

    // A preliminary check on Ye for nse solver.
    // This makes sure that input Ye is valid for this network.
    // If not valid, Ye is updated.

    amrex::Real ye_lo = 1.0_rt;
    amrex::Real ye_hi = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
        ye_lo = amrex::min(zion[n] * aion_inv[n], ye_lo);
        ye_hi = amrex::max(zion[n] * aion_inv[n], ye_hi);
    }

    if (state.y_e < ye_lo || state.y_e > ye_hi) {
#ifndef AMREX_USE_GPU
        if (state.y_e != -1.0_rt) {
            std::cout << "Warning: input Ye is " << state.y_e
                      << ", which is out of range. Ye is recalulated based on Composition"
                      << std::endl;
        }
#endif
        input_ye_is_valid = false;
    }

    if (!input_ye_is_valid) {
        // ensure Ye is valid
        // We assume the input state has rhoX

        state.y_e = 0.0_rt;
        for (int n = 0; n < NumSpec; ++n) {
            state.y_e += state.y[SFS+n] * zion[n] * aion_inv[n];
        }
        state.y_e /= state.rho;
    }

    if (input == nse_input_rty) {
        // Using (rho, T, Ye) as input

        return get_nse_state_rty(state, eps);

    } else if (input == nse_input_rey) {
        // Using (rho, e, Ye) as input

        return get_nse_state_rey(state, eps);

    } else {
        // When getting an unknown input

        amrex::Abort("Unknown NSE input mode");
        return state;
    }

}
#endif
