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
#include <actual_rhs.H>
#include <eos_composition.H>
#include <microphysics_sort.H>
#include <hybrj.H>
#include <cctype>
#include <algorithm>

using namespace nse_rp;

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
  amrex::Real dpf_dT;
  amrex::Real spin = 1.0_rt;

  // if we are doing drive_initial_convection, we want to use
  // the temperature that comes in through T_fixed

  Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

#ifndef NEW_NETWORK_IMPLEMENTATION
  auto tfactors = evaluate_tfactors(T_in);
#endif

  for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
      if (n == NSE_INDEX::h1_index) {
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


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_nse_exponent(T& nse_state) {
  // This function calculates a portion of the nse_state that dependents
  // on the chemical potential.
  // The input nse_state should be the result from get_nonexponent_nse_state

  // if we are doing drive_initial_convection, we want to use
  // the temperature that comes in through T_fixed

  Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
  Real exponent;

  // if we use a non-null screening method.
  // Get the required terms to calculate coulomb correction term, u_c

#if SCREEN_METHOD != 3

  // three unit-less constants for calculating coulomb correction term
  // see Calder 2007, doi:10.1086/510709 for more detail

  const amrex::Real A1 = -0.9052_rt;
  const amrex::Real A2 = 0.6322_rt;
  
  // A3 = -0.5 * sqrt(3) - A1/sqrt(A2)
  const amrex::Real A3 = 0.272433346667;

  // Find n_e for original state;
  const amrex::Real n_e = nse_state.rho * nse_state.y_e / C::m_u;
  const amrex::Real Gamma_e = C::q_e * C::q_e * std::cbrt(4.0_rt * M_PI * n_e / 3.0_rt)
                       / (C::k_B * T_in);
  amrex::Real gamma;
#endif

  amrex::Real u_c = 0.0_rt;
  
  for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
      if (n == NSE_INDEX::h1_index) {
          nse_state.xn[n] = 0.0;
          continue;
      }
#endif
      // term for calculating u_c

      // if use non-null screening, calculate the coulomb correction term
#if SCREEN_METHOD != 3
      gamma = std::pow(zion[n], 5.0_rt/3.0_rt) * Gamma_e;

      // chemical potential for coulomb correction

      u_c = C::k_B * T_in / C::Legacy::MeV2erg *
          (A1 * (std::sqrt(gamma * (A2 + gamma)) -
                 A2 * std::log(std::sqrt(gamma / A2) +
                               std::sqrt(1.0_rt + gamma / A2))) +
           2.0_rt * A3 * (std::sqrt(gamma) - std::atan(std::sqrt(gamma))));
#endif

      // find nse mass frac

      // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
      // iteration will make it happy again

      exponent = amrex::min(500.0_rt,
			    (zion[n] * nse_state.mu_p + (aion[n] - zion[n]) *
			     nse_state.mu_n - u_c + network::bion(n+1)) /
			    C::k_B / T_in * C::Legacy::MeV2erg);

      nse_state.xn[n] *= std::exp(exponent);
  }

  nse_state.y_e = 0.0_rt;

  for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
    if (n == NSE_INDEX::h1_index) {
      continue;
    }
#endif
    // constraint equation 1, mass fraction sum to 1

    nse_state.y_e += nse_state.xn[n] * zion[n] * aion_inv[n];
  }
}

// constraint equation

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fcn(Array1D<Real, 1, 2>& x, Array1D<Real, 1, 2>& fvec,
                const T& state, int& iflag) {
    // here state is the nse_state from get_nonexponent_nse_state
  
    T nse_state = state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    // Apply exponent component for calculating nse mass fractions

    apply_nse_exponent(nse_state);

    fvec(1) = -1.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::h1_index) {
            continue;
        }
#endif
        // constraint equation 1, mass fraction sum to 1

        fvec(1) += nse_state.xn[n];
    }

    // constraint equation 2, electron fraction should be the same

    fvec(2) = nse_state.y_e - state.y_e;

}

// constraint jacobian

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jcn(Array1D<Real, 1, 2>& x, Array2D<Real, 1, 2, 1, 2>& fjac,
                const T& state, int& iflag) {
    // here state is the nse_state from get_nonexponent_nse_state
  
    T nse_state = state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    apply_nse_exponent(nse_state);

    // evaluate jacobian of the constraint

    fjac(1, 1) = 0.0_rt;
    fjac(1, 2) = 0.0_rt;
    fjac(2, 1) = 0.0_rt;
    fjac(2, 2) = 0.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::h1_index) {
            continue;
        }
#endif
        fjac(1, 1) += nse_state.xn[n] * zion[n] / C::k_B / T_in * C::Legacy::MeV2erg ;
        fjac(1, 2) += nse_state.xn[n] * (aion[n] - zion[n]) / C::k_B / T_in * C::Legacy::MeV2erg;
        fjac(2, 1) += nse_state.xn[n] * zion[n] * zion[n] * aion_inv[n] / C::k_B / T_in * C::Legacy::MeV2erg;
        fjac(2, 2) += nse_state.xn[n] * zion[n] * (aion[n] - zion[n]) * aion_inv[n] / C::k_B / T_in * C::Legacy::MeV2erg;
    }

}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_hybrid_solver(T& state, amrex::Real eps=1.0e-10_rt) {
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

    outer_x(1) = state.mu_p;
    outer_x(2) = state.mu_n;

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
	hybrj(hj, state);

	fcn(hj.x, f, state, flag);

	if (std::abs(f(1)) < eps && std::abs(f(2)) < eps) {

	  state.mu_p = hj.x(1);
	  state.mu_n = hj.x(2);
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
    std::cout << "Temperature: " << state.T << std::endl;
    std::cout << "Density: " << state.rho << std::endl;
    std::cout << "Ye: " << state.y_e << std::endl;
    std::cout << "Initial mu_p and mu_n: " << state.mu_p << ", " << state.mu_n << std::endl;
#endif

    amrex::Error("failed to solve");
}

// A newton-raphson solver for finding nse state used for calibrating
// chemical potential of proton and neutron

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_nr_solver(T& state, amrex::Real eps=1.0e-10_rt) {
  // state is the nse_state from get_nonexponent_nse_state

  // whether nse solver converged or not

  bool converged = false;

  amrex::Array1D<amrex::Real, 1, 2> f;
  amrex::Array2D<amrex::Real, 1, 2, 1, 2> jac;
  amrex::Array1D<amrex::Real, 1, 2> x;
  int flag = 0;

  x(1) = state.mu_p;
  x(2) = state.mu_n;

  jcn(x, jac, state, flag);
  fcn(x, f, state, flag);

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
      state.mu_p = x(1);
      state.mu_n = x(2);
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

    jcn(x, jac, state, flag);
    fcn(x, f, state, flag);
  }

  if (!converged) {
    amrex::Error("NSE solver failed to converge!");
  }
}

// Get the NSE state;
template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_actual_nse_state(T& state, amrex::Real eps=1.0e-10_rt,
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
      std::cout << "Warning: input Ye is " << state.y_e << ", which is out of range. Ye is recalulated based on Composition" << std::endl;
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

  // Get nse_state without the exponent term

  auto nse_state = get_nonexponent_nse_state(state);

  // invoke newton-raphson or hybrj to solve chemical potential of proton and neutron
  // which are the exponent part of the nse mass fraction calculation

  if (use_hybrid_solver) {
    nse_hybrid_solver(nse_state, eps);
  }
  else {
    bool singular_network = true;
    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
      if (n == NSE_INDEX::h1_index) {
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

    nse_nr_solver(state, eps);
  }

  // Apply exponent for calculating nse mass fractions

  apply_nse_exponent(nse_state);

  // update mu_n and mu_p to input state

  state.mu_p = nse_state.mu_p;
  state.mu_n = nse_state.mu_n;
 
  return nse_state;
}
#endif
