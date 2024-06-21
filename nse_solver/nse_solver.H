#ifndef NSE_SOLVER_H
#define NSE_SOLVER_H

#include <fundamental_constants.H>
#include <AMReX_REAL.H>
#include <eos_type.H>
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
#include <approx_math.H>
#include <hybrj.H>
#include <screen.H>
#include <cctype>
#include <algorithm>


using namespace nse_rp;

template <typename T>
struct nse_solver_data
{
    T state;
    amrex::Array1D<amrex::Real, 1, NumSpec> u_c;
};


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_nonexponent_nse_state(const T& state) {
    // This function calculates a portion of the nse_state that is not dependent
    // on the chemical potential. This way we don't have to calculate the same
    // component repeatedly.

    // This function returns an incomplete nse_state where the xn lacks the exponent term.

    T nse_state = state;

    // set partition function and spin

    amrex::Real pf = 1.0_rt;
    [[maybe_unused]] amrex::Real dpf_dT;
    amrex::Real spin = 1.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

#ifndef NEW_NETWORK_IMPLEMENTATION
    auto tfactors = evaluate_tfactors(T_in);
#endif

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            nse_state.xn[n] = 0.0;
            continue;
        }
#endif

        // fill partition function and get spin

#ifndef NEW_NETWORK_IMPLEMENTATION
        spin = get_spin_state(n+1);
        get_partition_function(n+1, tfactors, pf, dpf_dT);
#endif

        // find nse mass frac without the exponent term.

        Real power = 2.0 * M_PI * network::mion(n+1) *
            C::k_B * T_in / (C::hplanck * C::hplanck);

        nse_state.xn[n] = network::mion(n+1) * pf * spin / state.rho *
            std::sqrt(power * power * power);
    }

    return nse_state;
}

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void compute_coulomb_contribution(amrex::Array1D<amrex::Real, 1, NumSpec>& u_c,
                                const T& state) {
    // This function computes the coulomb contribution 1D array

    // if we use chabrier1998 screening
    // Get the required terms to calculate coulomb correction term, u_c

    amrex::Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    //
    // Find n_e for original state;
    // Note that y_e depends on the mass fraction,
    // but we use the coulomb correction to compute the mass fraction
    // So here y_e is simply the actual y_e we want to achieve.
    // so we just treat u_c as a constant.
    //

    const amrex::Real n_e = state.rho * state.y_e / C::m_u;
    const amrex::Real Gamma_e = C::q_e * C::q_e *
        std::cbrt(4.0_rt * M_PI * n_e / 3.0_rt) / (C::k_B * T_in);
    amrex::Real gamma;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // term for calculating u_c

        gamma = std::pow(zion[n], 5.0_rt/3.0_rt) * Gamma_e;

        // chemical potential for coulomb correction
        // see appendix of Calder 2007, doi:10.1086/510709 for more detail

        // reuse existing implementation from screening routine
        amrex::Real f, df;
        constexpr int do_T_derivatives = 0;
        chabrier1998_helmholtz_F<do_T_derivatives>(gamma, 0.0_rt, f, df);

        //
        // Here u_c is a dimensionless quantity.
        // Otherwise:
        // u_c = C::k_B * T_in / C::Legacy::MeV2erg * f;
        //

        u_c(n+1) = f;
    }
}
#endif

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_nse_exponent(T& nse_state,
                        const amrex::Array1D<amrex::Real, 1, NumSpec>& u_c) {
    // This function applies the nse exponent which depends on
    // the chemical potential.

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    amrex::Real exponent;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            nse_state.xn[n] = 0.0;
            continue;
        }
#endif

        // find nse mass frac

        // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
        // iteration will make it happy again

        exponent = amrex::min(500.0_rt,
                              (zion[n] * nse_state.mu_p + (aion[n] - zion[n]) *
                               nse_state.mu_n + network::bion(n+1)) /
                              C::k_B / T_in * C::Legacy::MeV2erg - u_c(n+1));

        nse_state.xn[n] *= fast_exp(exponent);
    }

    nse_state.y_e = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // constraint equation 1, mass fraction sum to 1

        nse_state.y_e += nse_state.xn[n] * zion[n] * aion_inv[n];
    }
}

// constraint equation

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fcn(Array1D<Real, 1, 2>& x, Array1D<Real, 1, 2>& fvec,
         const nse_solver_data<T>& state_data, int& iflag) {
    // here state is the nse_state from get_nonexponent_nse_state

    amrex::ignore_unused(iflag);

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    // Apply exponent component for calculating nse mass fractions

    apply_nse_exponent(nse_state, state_data.u_c);

    fvec(1) = -1.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // constraint equation 1, mass fraction sum to 1

        fvec(1) += nse_state.xn[n];
    }

    // constraint equation 2, electron fraction should be the same

    fvec(2) = nse_state.y_e - state_data.state.y_e;

}

// constraint jacobian

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jcn(Array1D<Real, 1, 2>& x, Array2D<Real, 1, 2, 1, 2>& fjac,
         const nse_solver_data<T>& state_data, int& iflag) {
    // here state is the nse_state from get_nonexponent_nse_state

    amrex::ignore_unused(iflag);

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    apply_nse_exponent(nse_state, state_data.u_c);

    // evaluate jacobian of the constraint

    fjac(1, 1) = 0.0_rt;
    fjac(1, 2) = 0.0_rt;
    fjac(2, 1) = 0.0_rt;
    fjac(2, 2) = 0.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        fjac(1, 1) += nse_state.xn[n] * zion[n] / C::k_B / T_in * C::Legacy::MeV2erg ;
        fjac(1, 2) += nse_state.xn[n] * (aion[n] - zion[n]) / C::k_B / T_in * C::Legacy::MeV2erg;
        fjac(2, 1) += nse_state.xn[n] * zion[n] * zion[n] * aion_inv[n] / C::k_B / T_in * C::Legacy::MeV2erg;
        fjac(2, 2) += nse_state.xn[n] * zion[n] * (aion[n] - zion[n]) * aion_inv[n] / C::k_B / T_in * C::Legacy::MeV2erg;
    }

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_hybrid_solver(nse_solver_data<T>& state_data,
                       amrex::Real eps=1.0e-10_rt) {
    // state is the nse_state from get_nonexponent_nse_state

    hybrj_t<2> hj;

    // we'll take x[1] = mu_p, x[2] = mu_n

    hj.xtol = eps;
    hj.mode = 1;

    // random flag number for evaluation in for loop;

    int flag = 0;

    // Fine-tune variables

    amrex::Real dx;
    bool is_pos_new;
    bool is_pos_old = false;

    amrex::Array1D<amrex::Real, 1, 2> f;
    amrex::Array1D<amrex::Real, 1, 2> outer_x;
    amrex::Array1D<amrex::Real, 1, 2> inner_x;

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

            // hybrj<2, T>(hj, state, fcn_hybrid<T>, jcn_hybrid<T>);
            hybrj(hj, state_data);

            fcn(hj.x, f, state_data, flag);

            if (std::abs(f(1)) < eps && std::abs(f(2)) < eps) {

                state_data.state.mu_p = hj.x(1);
                state_data.state.mu_n = hj.x(2);
                return;
            }

            if (f(1) > 0.0_rt && f(2) > 0.0_rt) {
                is_pos_new = true;
            }
            else{
                is_pos_new = false;
            }

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

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_nr_solver(nse_solver_data<T>& state_data,
                   amrex::Real eps=1.0e-10_rt) {
    // state_data is the state_data after from
    // get_nonexponent_nse_state and compute_coulomb_contribution

    // whether nse solver converged or not

    bool converged = false;

    amrex::Array1D<amrex::Real, 1, 2> f;
    amrex::Array2D<amrex::Real, 1, 2, 1, 2> jac;
    amrex::Array1D<amrex::Real, 1, 2> x;
    int flag = 0;

    x(1) = state_data.state.mu_p;
    x(2) = state_data.state.mu_n;

    jcn(x, jac, state_data, flag);
    fcn(x, f, state_data, flag);

    // store determinant for finding inverse jac
    amrex::Real det;

    // store inverse jacobian
    amrex::Array2D<amrex::Real, 1, 2, 1, 2> inverse_jac;

    // difference in chemical potential of proton and neutron
    amrex::Real d_mu_p = std::numeric_limits<Real>::max();
    amrex::Real d_mu_n = std::numeric_limits<Real>::max();

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

        inverse_jac(1, 1) = jac(2, 2) / scale_fac / det;
        inverse_jac(1, 2) = -jac(1, 2) / scale_fac / det;
        inverse_jac(2, 1) = -jac(2, 1) / scale_fac / det;
        inverse_jac(2, 2) = jac(1, 1) / scale_fac / det;

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

        jcn(x, jac, state_data, flag);
        fcn(x, f, state_data, flag);
    }

    if (!converged) {
        amrex::Error("NSE solver failed to converge!");
    }
}

// Get the NSE state;
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_actual_nse_state(T& state, amrex::Real eps=1.0e-6_rt,
                       bool input_ye_is_valid=false) {

    // Check whether input ye is actually valid
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
#ifdef STRANG
        composition(state);
#else
        state.y_e = 0.0_rt;

        for (int n = 0; n < NumSpec; ++n) {
            state.y_e += state.y[SFS+n] * zion[n] * aion_inv[n] / state.rho;
        }
#endif
    }

    nse_solver_data<T> state_data = {state, {0.0_rt}};

    // Get nse_state without the exponent term

    state_data.state = get_nonexponent_nse_state(state);

    // if use chabrier1998 screening, calculate the coulomb correction term

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    compute_coulomb_contribution(state_data.u_c, state);
#endif

    // invoke newton-raphson or hybrj to solve chemical potential of proton and neutron
    // which are the exponent part of the nse mass fraction calculation

    if (use_hybrid_solver) {
        nse_hybrid_solver(state_data, eps);
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

        if (singular_network == true) {
            amrex::Error("This network always results in singular jacobian matrix, thus can't find nse mass fraction using nr!");
        }

        nse_nr_solver(state_data, eps);
    }

    // Apply exponent for calculating nse mass fractions

    apply_nse_exponent(state_data.state, state_data.u_c);

    // update mu_n and mu_p to input state

    state.mu_p = state_data.state.mu_p;
    state.mu_n = state_data.state.mu_n;

    return state_data.state;
}
#endif
