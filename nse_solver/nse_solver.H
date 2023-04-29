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
#include <eos_composition.H>
#include <microphysics_sort.H>
#include <nse_index.H>
#include <hybrj.H>
#include <cctype>
#include <algorithm>

AMREX_INLINE
void find_rate_indices(amrex::Array1D<int, 1, 3>& reactant_indices,
		       amrex::Array1D<int, 1, 3>& product_indices, const std::string& rate_name){

  // fill the corresponding reactant and product indices of a given rate name
  // used during initialization

  for (int n = 1; n <= 3; ++n){
    reactant_indices(n) = -1;
    product_indices(n) = -1;
  }

  std::string rate_name_mod = rate_name;

  // Process rate_name to find the reactants and products

  // replace p with h1 to match with short_spec_name_cxx

  size_t pos;

  while (true){
    pos = rate_name_mod.find("p_");

    if (pos == std::string::npos){
      break;
    }

    rate_name_mod.replace(pos, 2, "h1_");
  }

  // Separate product name from reactant names

  pos = rate_name_mod.find("__");
  std::string reactant_names = rate_name_mod.substr(0, pos);
  std::string product_names = rate_name_mod.substr(pos+2);

  // Remove suffix in the product name

  std::string suffix;

  pos = product_names.find("__");
  if (pos != std::string::npos){
    suffix = product_names.substr(pos);
    product_names.erase(pos);
  }

  // check whether if this rate is removed from the network
  // let its reactant and product all to be -1 if rate is removed

  pos = suffix.find("removed");
  if (pos != std::string::npos){
    return;
  }

  std::string _temp;

  // Find the index of reactants and products in short_spec_names

  int i = 1;
  while ((pos = reactant_names.find('_')) != std::string::npos){

    _temp = reactant_names.substr(0, pos);
    // aprox network need capitalize
    // but we only work with pynucastro networks for now.

    // cap_first_letter(reactant);

    for (int n = 0; n < NumSpec; ++n){
      if (short_spec_names_cxx[n] == _temp){
	reactant_indices(i) = n;
	++i;
      }
    }

    reactant_names.erase(0, pos+1);
  }

  // cap_first_letter(reactant_names);
  for (int n = 0; n < NumSpec; ++n){
    if (short_spec_names_cxx[n] == reactant_names){
      reactant_indices(i) = n;
    }
  }


  i = 1;
  while ((pos = product_names.find('_')) != std::string::npos){

    _temp = product_names.substr(0, pos);
    // cap_first_letter(product);

    for (int n = 0; n < NumSpec; ++n){
      if (short_spec_names_cxx[n] == _temp){
	product_indices(i) = n;
	++i;
      }
    }

    product_names.erase(0, pos+1);
  }

  // cap_first_letter(product_names);

  for (int n = 0; n < NumSpec; ++n){
    if (short_spec_names_cxx[n] == product_names){
      product_indices(i) = n;
    }
  }

  // Sort rate_indices.

  sort_Array1D<1, 3>(reactant_indices);
  sort_Array1D<1, 3>(product_indices);

}


AMREX_INLINE
void init_nse_net() {
  // Initialization for nse, check whether network is valid and stores indices

  for (int n = 0; n < NumSpec; ++n){

    // store index of photoionization proton by looking for Z=1 but
    // exclude h1 which is always at index=0 for aprox networks
#ifdef NEW_NETWORK_IMPLEMENTATION
    if (zion[n] == 1 && aion[n] == 1){
      if (n == 0){
	NSE_INDEX::h1_index = n;
      }
      else {
	NSE_INDEX::p_index = n;
      }
    }
#else
    if (zion[n] == 1 && aion[n]== 1){
      NSE_INDEX::p_index = n;
    }
#endif
    else if (zion[n] == 0){
      NSE_INDEX::n_index = n;
    }
    else if (zion[n] == 2 && aion[n] == 4){
      NSE_INDEX::he4_index = n;
    }
  }

  // Check if network results in singular jacobian first, require at
  // least one nuclei that nuclei.Z != nuclei.N.  Some examples include
  // aprox13 and iso7. If use hybrj solver, this is no longer an issue.
  if (!use_hybrid_solver){
    bool singular_network = true;
    for (int n = 0; n < NumSpec; ++n){
      if (n == NSE_INDEX::h1_index){
	continue;
      }
      if (zion[n] != aion[n] - zion[n]){
	singular_network = false;
      }
    }

    if (singular_network == true){
      amrex::Error("This network always results in singular jacobian matrix, thus can't find nse mass fraction using nr!");
    }
  }

#ifndef NEW_NETWORK_IMPLEMENTATION
  // rate names initialization for non aprox networks...

  // Fill in rate_indices<int, 1, NumRates, 1, 6>
  // used for referencing rate by nuclei indices

  amrex::Array1D<int, 1, 3> reactant_indices;
  amrex::Array1D<int, 1, 3> product_indices;
  std::string current_rate_name;

  for (int n = 1; n <= Rates::NumRates; ++n){

    // get current name of the rate
    current_rate_name = Rates::rate_names[n];

    // fill in reactant_indices and product indices
    find_rate_indices(reactant_indices, product_indices, current_rate_name);

    // fill in rate_indices
    for (int i = 1; i <= 3; ++i){
      NSE_INDEX::rate_indices(n, i) = reactant_indices(i);
      NSE_INDEX::rate_indices(n, i+3) = product_indices(i);
    }

  }
#endif

  // make sure nse_molar_independent is used if nse_skip_molar is enabled

  if (nse_skip_molar) {
    nse_molar_independent = true;
  }
  
  NSE_INDEX::initialized = true;
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_nse_state(const T& state)
{
  // This function finds the nse state given the burn state or eos state

  // three unit-less constants for calculating coulomb correction term
  // See Calder 2007, doi:10.1086/510709 paper for more detail


  // Store nse_state
  T nse_state;

  // Get the required terms to calculate coulomb correction term, u_c
  // if we use a non-null screening method.
  
#if SCREEN_METHOD != 3
  const amrex::Real A1 = -0.9052_rt;
  const amrex::Real A2 = 0.6322_rt;
  const amrex::Real A3 = -0.5_rt * std::sqrt(3.0_rt) - A1 / std::sqrt(A2);

  // Find n_e for original state;
  const amrex::Real n_e = state.rho * state.y_e / C::m_u;
  amrex::Real gamma;
#endif

  amrex::Real u_c = 0.0_rt;
  
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

  for (int n = 0; n < NumSpec; ++n){
      if (n == NSE_INDEX::h1_index){
          nse_state.xn[n] = 0.0;
          continue;
      }

      // term for calculating u_c

      // if use non-null screening, calculate the coulomb correction term

#if SCREEN_METHOD != 3
      gamma = std::pow(zion[n], 5.0_rt/3.0_rt) *
          C::q_e * C::q_e * std::cbrt(4.0_rt * M_PI * n_e / 3.0_rt) /
          (C::k_B * T_in);

      // chemical potential for coulomb correction

      u_c = C::k_B * T_in / C::Legacy::MeV2erg *
          (A1 * (std::sqrt(gamma * (A2 + gamma)) -
                 A2 * std::log(std::sqrt(gamma / A2) +
                               std::sqrt(1.0_rt + gamma / A2))) +
           2.0_rt * A3 * (std::sqrt(gamma) - std::atan(std::sqrt(gamma))));
#endif

      // fill partition function and get spin

#ifndef NEW_NETWORK_IMPLEMENTATION
      spin = get_spin_state(n+1);
      get_partition_function(n+1, tfactors, pf, dpf_dT);
#endif

      // find nse mass frac

      // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
      // iteration will make it happy again

      Real exponent = amrex::min(500.0_rt,
                                 (zion[n] * state.mu_p + (aion[n] - zion[n]) * state.mu_n
                                  - u_c + network::bion(n+1)) /
                                 C::k_B / T_in * C::Legacy::MeV2erg);

      nse_state.xn[n] = network::mion(n+1) * pf * spin / state.rho *
          std::pow(2.0 * M_PI * network::mion(n+1) *
                   C::k_B * T_in / std::pow(C::hplanck, 2.0_rt), 3.0_rt/2.0_rt) *
          std::exp(exponent);
  }

  nse_state.y_e = 0.0_rt;

  for (int n = 0; n < NumSpec; ++n){
    if (n == NSE_INDEX::h1_index){
      continue;
    }

    // constraint equation 1, mass fraction sum to 1

    nse_state.y_e += nse_state.xn[n] * zion[n] * aion_inv[n];
  }

  nse_state.T = T_in;
  nse_state.rho = state.rho;
  nse_state.mu_p = state.mu_p;
  nse_state.mu_n = state.mu_n;
  
  return nse_state;
}



// constraint equation

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fcn_hybrid(Array1D<Real, 1, 2>& x, Array1D<Real, 1, 2>& fvec,
                const T& state, int& iflag) {

    T current_state = state;
    current_state.mu_p = x(1);
    current_state.mu_n = x(2);

    auto nse_state = get_nse_state(current_state);

    fvec(1) = -1.0_rt;

    for (int n = 0; n < NumSpec; ++n){
        if (n == NSE_INDEX::h1_index){
            continue;
        }

        // constraint equation 1, mass fraction sum to 1

        fvec(1) += nse_state.xn[n];
    }

    // constraint equation 2, electron fraction should be the same

    fvec(2) = nse_state.y_e - state.y_e;

}

// constraint jacobian

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jcn_hybrid(Array1D<Real, 1, 2>& x, Array2D<Real, 1, 2, 1, 2>& fjac,
                const T& state, int& iflag) {

    T current_state = state;
    current_state.mu_p = x(1);
    current_state.mu_n = x(2);

    auto nse_state = get_nse_state(current_state);

    // evaluate jacobian of the constraint

    fjac(1, 1) = 0.0_rt;
    fjac(1, 2) = 0.0_rt;
    fjac(2, 1) = 0.0_rt;
    fjac(2, 2) = 0.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    for (int n = 0; n < NumSpec; ++n){
        if (n == NSE_INDEX::h1_index){
            continue;
        }

        fjac(1, 1) += nse_state.xn[n] * zion[n] / C::k_B / T_in * C::Legacy::MeV2erg ;
        fjac(1, 2) += nse_state.xn[n] * (aion[n] - zion[n]) / C::k_B / T_in * C::Legacy::MeV2erg;
        fjac(2, 1) += nse_state.xn[n] * zion[n] * zion[n] * aion_inv[n] / C::k_B / T_in * C::Legacy::MeV2erg;
        fjac(2, 2) += nse_state.xn[n] * zion[n] * (aion[n] - zion[n]) * aion_inv[n] / C::k_B / T_in * C::Legacy::MeV2erg;
    }

}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_hybrid_solver(T& state, amrex::Real eps=1.0e-10_rt) {

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

    for (int i = 0; i < 20; ++i){

      dx = 0.5_rt;
      inner_x(1) = outer_x(1);
      inner_x(2) = outer_x(2);

      for (int j = 0; j < 20; ++j){

	hj.x(1) = inner_x(1);
	hj.x(2) = inner_x(2);

	hybrj(hj, state, fcn_hybrid<T>, jcn_hybrid<T>);

	fcn_hybrid(hj.x, f, state, flag);

	if (std::abs(f(1)) < eps && std::abs(f(2)) < eps){

	  state.mu_p = hj.x(1);
	  state.mu_n = hj.x(2);
	  return;
	}

	if (f(1) > 0.0_rt && f(2) > 0.0_rt){
	  is_pos_new = true;
	}
	else{
	  is_pos_new = false;
	}

	if (is_pos_old != is_pos_new){
	  dx *= 0.8_rt;
	}

	if (is_pos_new){
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

    amrex::Error("failed to solve");
}

// A newton-raphson solver for finding nse state used for calibrating
// chemical potential of proton and neutron

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_nr_solver(T& state, amrex::Real eps=1.0e-10_rt) {

  bool converged = false;                                     // whether nse solver converged or not

  amrex::Array1D<amrex::Real, 1, 2> f;
  amrex::Array2D<amrex::Real, 1, 2, 1, 2> jac;
  amrex::Array1D<amrex::Real, 1, 2> x;
  int flag = 0;

  x(1) = state.mu_p;
  x(2) = state.mu_n;

  jcn_hybrid(x, jac, state, flag);
  fcn_hybrid(x, f, state,flag);

  amrex::Real det;                                            // store determinant for finding inverse jac
  amrex::Array2D<amrex::Real, 1, 2, 1, 2> inverse_jac;        // store inverse jacobian
  amrex::Real d_mu_p = std::numeric_limits<Real>::max();      // difference in chemical potential of proton
  amrex::Real d_mu_n = std::numeric_limits<Real>::max();      // difference in chemical potential of neutron

  // begin newton-raphson
  for (int i = 0; i < max_nse_iters; ++i){

    // check if current state fulfills constraint equation

    if (std::abs(d_mu_p) < eps * std::abs(x(1)) &&
        std::abs(d_mu_n) < eps * std::abs(x(2))){
      converged = true;
      state.mu_p = x(1);
      state.mu_n = x(2);
      break;
    }

    // Find the max of the jacobian used for scaling determinant to prevent digit overflow

    auto scale_fac = amrex::max(jac(2,2), amrex::max(jac(2,1), amrex::max(jac(1,1), jac(1,2))));

    // if jacobians are small, then no need for scaling

    if (scale_fac < 1.0e150){
      scale_fac = 1.0_rt;
    }

    // Specific inverse 2x2 matrix, perhaps can write a function for solving n systems of equations.

    det = jac(1, 1) / scale_fac * jac(2, 2) - jac(1, 2) / scale_fac * jac(2, 1);

    // check if determinant is 0

    if (det == 0.0_rt){
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

    if (std::abs(d_mu_p) > 1.0e3_rt || std::abs(d_mu_n) > 1.0e3_rt){
      amrex::Error("Not making good progress, breaking");
    }

    // update new solution

    x(1) += d_mu_p;
    x(2) += d_mu_n;

    // check whether solution results in nan

    if (std::isnan(x(1)) || std::isnan(x(2))){
	amrex::Error("Nan encountered, likely due to overflow in digits or not making good progress");
      }

    // update constraint

    jcn_hybrid(x, jac, state, flag);
    fcn_hybrid(x, f, state, flag);
  }

  if (!converged){
    amrex::Error("NSE solver failed to converge!");
  }
}

// Get the NSE state;
template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_actual_nse_state(T& state, amrex::Real eps=1.0e-10_rt,
                       bool input_ye_is_valid=false){

  // Check whether initialized or not
  if (!NSE_INDEX::initialized){
    amrex::Error("NSE uninitialized! Need to call init_nse_net() ");
  }
  // Check whether input ye is actually valid
  amrex::Real ye_lo = 1.0_rt;
  amrex::Real ye_hi = 0.0_rt;

  for (int n = 0; n < NumSpec; ++n){
    ye_lo = amrex::min(zion[n] * aion_inv[n], ye_lo);
    ye_hi = amrex::max(zion[n] * aion_inv[n], ye_hi);
  }

  if (state.y_e < ye_lo || state.y_e > ye_hi){
    input_ye_is_valid = false;
  }

  if (!input_ye_is_valid) {
    // ensure Ye is valid
    composition(state);
  }

  // invoke newton-raphson to solve chemical potential of proton and neutron
  if (use_hybrid_solver) {
    nse_hybrid_solver(state, eps);
  }
  else{
    nse_nr_solver(state, eps);
  }

  auto nse_state = get_nse_state(state);

  return nse_state;
}
#endif
