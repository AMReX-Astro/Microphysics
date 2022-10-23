#ifndef _actual_nse_H_
#define _actual_nse_H_

#include <fundamental_constants.H>
#include <AMReX_REAL.H>
#include <eos_type.H>
#include <network.H>
#include <burner.H>
#include <extern_parameters.H>
#include <cmath>
#include <AMReX_Array.H>
#include <actual_network.H>
#include <eos_composition.H>
#include <nse_index.H>
#include <hybrj.H>
#include <cctype>
#include <algorithm>


template <int l, int m>
void sort_Array1D(amrex::Array1D<int, l, m>& array_1d){

  // a simple bubble sort, for sorting Array1D
  // assumes index start from int l, and end at int m

  int scratch;
  
  for (int i = l; i <= m; ++i){
    for (int j = l; j <= m-i; ++j){
      if (array_1d(j) > array_1d(j+1)){
	scratch = array_1d(j);
	array_1d(j) = array_1d(j+1);
	array_1d(j+1) = scratch;
      }
    }
  }
}


template <int l, int m>
bool in_single_group(const amrex::Array1D<int, l, m>& group_ind){

  // This function checks whether a 1d array has the same element overall
  
  for (int n = l; n <= m; ++n){
    if (group_ind(l) != group_ind(n)){
      return false;
    }
  }

  return true;
}



void find_rate_indices(amrex::Array1D<int, 1, 3>& reactant_indices,
		       amrex::Array1D<int, 1, 3>& product_indices, const std::string& rate_name){

  // fill the corresponding reactant and product indices of a given rate name
  
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
  
  pos = product_names.find("__");
  if (pos != std::string::npos){
    std::string suffix = product_names.substr(pos);
    product_names.erase(pos);
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


void init_actual_nse(){
  
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
    if (zion[n] ==1 && aion[n] == 1){
      NSE_INDEX::h1_index = n;
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
  // aprox13 and iso7

  bool singular_network = true;
  for (int n = 0; n < NumSpec; ++n){
    if (n == NSE_INDEX::p_index){
      continue;
    }
    if (zion[n] != aion[n] - zion[n]){
      singular_network = false;
    }
  }

  if (singular_network == true){
    amrex::Error("This network always results in singular jacobian matrix, thus can't find nse mass fraction!");
  }

  // rate names initialization...

  // Fill in rate_indices<int, 1, NumRates, 1, 7>
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
    
    NSE_INDEX::rate_indices(n, 7) = n;
  }

  
  NSE_INDEX::initialized = true;
}


// struct to store output from constraint equations and jacobian for
// newton-raphson
struct Newton_inputs{
  amrex::Array1D<amrex::Real, 0, 1> eqs;
  amrex::Array2D<amrex::Real, 0, 1, 0, 1> jac;
};

template <typename T>
T get_nse_state(const T& state)
{
  // This function finds the nse state given the burn state or eos state

  // three unit-less constants for calculating coulomb correction term
  // See Calder 2007, doi:10.1086/510709 paper for more detail

  const amrex::Real A1 = -0.9052_rt;
  const amrex::Real A2 = 0.6322_rt;
  const amrex::Real A3 = -0.5_rt * std::sqrt(3.0_rt) - A1 / std::sqrt(A2);

  // Store nse_state
  T nse_state;

  // Find n_e for original state;

  const amrex::Real n_e = state.rho * state.y_e / C::m_u;
  amrex::Real gamma;
  amrex::Real u_c;

  // Need partition function, set it to 1 for now.

  amrex::Real partition_function = 1.0_rt;

  for (int n = 0; n < NumSpec; ++n){
      if (n == NSE_INDEX::p_index){
          continue;
      }

      // term for calculating u_c

      gamma = std::pow(zion[n], 5.0_rt/3.0_rt) *
          C::q_e * C::q_e * std::cbrt(4.0_rt * M_PI * n_e / 3.0_rt) /
          (C::k_B * state.T);

      // chemical potential for coulomb correction

      u_c = C::k_B * state.T / C::Legacy::MeV2erg *
          (A1 * (std::sqrt(gamma * (A2 + gamma)) -
                 A2 * std::log(std::sqrt(gamma / A2) +
                               std::sqrt(1.0_rt + gamma / A2))) +
           2.0_rt * A3 * (std::sqrt(gamma) - std::atan(std::sqrt(gamma))));


      // find nse mass frac

      // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
      // iteration will make it happy again

      Real exponent = amrex::min(500.0_rt,
                                 (zion[n] * state.mu_p + (aion[n] - zion[n]) * state.mu_n
                                  - u_c + network::bion(n+1)) /
                                 C::k_B / state.T * C::Legacy::MeV2erg);

      nse_state.xn[n] = network::mion(n+1) * partition_function / state.rho *
          std::pow(2.0 * M_PI * network::mion(n+1) *
                   C::k_B * state.T / std::pow(C::hplanck, 2.0_rt), 3.0_rt/2.0_rt) *
          std::exp(exponent);
  }

  nse_state.T = state.T;
  nse_state.rho = state.rho;

  return nse_state;

}


template <typename T>
Newton_inputs nse_constraint(const T& state){
    // This functions finds the constraint equations and jacobian used
    // for calculating nse.

    Newton_inputs nse_inputs;

    // calculate the nse state based on initial conditions of mu_p and
    // mu_n which are chemical potential of proton and neutron

    auto nse_state = get_nse_state(state);

    // Now find constraint equations

    nse_state.y_e = 0.0_rt;

    nse_inputs.eqs(0) = -1.0_rt;

    for (int n = 0; n < NumSpec; ++n){
        if (n == NSE_INDEX::p_index){
            continue;
        }

        // constraint equation 1, mass fraction sum to 1

        nse_inputs.eqs(0) += nse_state.xn[n];
        nse_state.y_e += nse_state.xn[n] * zion[n] * aion_inv[n];
    }

    // constraint equation 2, electron fraction should be the same

    nse_inputs.eqs(1) = nse_state.y_e - state.y_e;

    // evaluate jacobian of the constraint

    nse_inputs.jac(0,0) = 0.0_rt;
    nse_inputs.jac(0,1) = 0.0_rt;
    nse_inputs.jac(1,0) = 0.0_rt;
    nse_inputs.jac(1,1) = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n){
        if (n == NSE_INDEX::p_index){
            continue;
        }

        nse_inputs.jac(0,0) += nse_state.xn[n] * zion[n] / C::k_B / state.T * C::Legacy::MeV2erg ;
        nse_inputs.jac(0,1) += nse_state.xn[n] * (aion[n] - zion[n]) / C::k_B / state.T * C::Legacy::MeV2erg;
        nse_inputs.jac(1,0) += nse_state.xn[n] * zion[n] * zion[n] * aion_inv[n] / C::k_B / state.T * C::Legacy::MeV2erg;
        nse_inputs.jac(1,1) += nse_state.xn[n] * zion[n] * (aion[n] - zion[n]) * aion_inv[n] / C::k_B / state.T * C::Legacy::MeV2erg;
    }

    return nse_inputs;
}


// for the hybrid Powell solver

template<typename T>
void fcn_hybrid(Array1D<Real, 1, 2>& x, Array1D<Real, 1, 2>& fvec,
                const T& state, int& iflag) {

    T current_state = state;
    current_state.mu_p = x(1);
    current_state.mu_n = x(2);

    auto nse_state = get_nse_state(current_state);

    nse_state.y_e = 0.0_rt;

    fvec(1) = -1.0_rt;

    for (int n = 0; n < NumSpec; ++n){
        if (n == NSE_INDEX::p_index){
            continue;
        }

        // constraint equation 1, mass fraction sum to 1

        fvec(1) += nse_state.xn[n];
        nse_state.y_e += nse_state.xn[n] * zion[n] * aion_inv[n];
    }

    // constraint equation 2, electron fraction should be the same

    fvec(2) =  nse_state.y_e - state.y_e;

}

template<typename T>
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

    for (int n = 0; n < NumSpec; ++n){
        if (n == NSE_INDEX::p_index){
            continue;
        }

        fjac(1, 1) += nse_state.xn[n] * zion[n] / C::k_B / state.T * C::Legacy::MeV2erg ;
        fjac(1, 2) += nse_state.xn[n] * (aion[n] - zion[n]) / C::k_B / state.T * C::Legacy::MeV2erg;
        fjac(2, 1) += nse_state.xn[n] * zion[n] * zion[n] * aion_inv[n] / C::k_B / state.T * C::Legacy::MeV2erg;
        fjac(2, 2) += nse_state.xn[n] * zion[n] * (aion[n] - zion[n]) * aion_inv[n] / C::k_B / state.T * C::Legacy::MeV2erg;
    }

}

template<typename T>
void nse_hybrid_solver(T& state, amrex::Real eps=1.0e-3_rt) {

    hybrj_t<2> hj;

    // we'll take x[1] = mu_p, x[2] = mu_n

    hj.x(1) = state.mu_p;
    hj.x(2) = state.mu_n;

    hj.xtol = eps;
    hj.mode = 2;

    for (int j = 1; j <= 2; ++j) {
        hj.diag(j) = 1.0_rt;
    }

    hybrj(hj, state, fcn_hybrid<T>, jcn_hybrid<T>);

    if (hj.info != 1) {
        amrex::Error("failed to solve");
    }

    state.mu_p = hj.x(1);
    state.mu_n = hj.x(2);

}

// A newton-raphson solver for finding nse state used for calibrating
// chemical potential of proton and neutron
template<typename T>
void nse_nr_solver(T& state, amrex::Real eps=1.0e-3_rt) {

  bool converged = false;                                     // whether nse solver converged or not

  Newton_inputs f = nse_constraint(state);                    // get constraint eqs and jacobian

  amrex::Real det;                                            // store determinant for finding inverse jac
  decltype(f.jac) inverse_jac;                                // store inverse jacobian
  amrex::Real d_mu_p = std::numeric_limits<Real>::max();                                // difference in chemical potential of proton
  amrex::Real d_mu_n = std::numeric_limits<Real>::max();                                // difference in chemical potential of neutron

  // begin newton-raphson
  for (int i = 0; i < max_nse_iters; ++i){

    // check if current state fulfills constraint equation
//    if (std::abs(f.eqs(0)) < eps && std::abs(f.eqs(1)) < eps){
    if (std::abs(d_mu_p) < eps * std::abs(state.mu_p) &&
        std::abs(d_mu_n) < eps * std::abs(state.mu_n)){
      converged = true;
      break;
    }

    // Find the max of the jacobian used for scaling determinant to prevent digit overflow
    auto scale_fac = amrex::max(f.jac(1,1),amrex::max(f.jac(1,0), amrex::max(f.jac(0,0), f.jac(0,1))));

    // if jacobians are small, then no need for scaling
    if (scale_fac < 1.0e150){
      scale_fac = 1.0_rt;
    }

    // Specific inverse 2x2 matrix, perhaps can write a function for solving n systems of equations.
    det = f.jac(0, 0) / scale_fac * f.jac(1, 1) - f.jac(0, 1) / scale_fac * f.jac(1, 0);

    // check if determinant is 0
    if (det == 0.0_rt){
      amrex::Error("Jacobian is a singular matrix! Try a different initial guess!");
    }

    // find inverse jacobian
    inverse_jac(0, 0) = f.jac(1,1) / scale_fac / det;
    inverse_jac(0, 1) = -f.jac(0,1) / scale_fac / det;
    inverse_jac(1, 0) = -f.jac(1,0) / scale_fac / det;
    inverse_jac(1, 1) = f.jac(0,0) / scale_fac / det;

    // find the difference
    d_mu_p = -(f.eqs(0) * inverse_jac(0,0) + f.eqs(1) * inverse_jac(0,1));
    d_mu_n = -(f.eqs(0) * inverse_jac(1,0) + f.eqs(1) * inverse_jac(1,1));

    // if diff goes beyond 1.0e3_rt, likely that its not making good progress..
    if (std::abs(d_mu_p) > 1.0e3_rt or std::abs(d_mu_n) > 1.0e3_rt){
      amrex::Error("Not making good progress, breaking");
    }

    // update new solution
    state.mu_p += d_mu_p;
    state.mu_n += d_mu_n;

    // check whether solution results in nan
    if (std::isnan(state.mu_p) or std::isnan(state.mu_n)){
      amrex::Error("Nan encountered, likely due to overflow in digits or not making good progress");
    }

    // update constraint
    f = nse_constraint(state);
  }

  if (!converged){
    amrex::Error("NSE solver failed to converge!");
  }
}

// Get the NSE state;
template<typename T>
T get_actual_nse_state(T& state, amrex::Real eps=1.0e-4_rt, bool input_ye_is_valid=false){

  // Check whether initialized or not
  if (!NSE_INDEX::initialized){
    amrex::Error("NSE uninitialized! Need to call init_actual_nse() ");
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
  } else {
      nse_nr_solver(state, eps);
  }

  // get the nse_state
  T nse_state = get_nse_state(state);
  return nse_state;
}


// First check to see if we're in the ballpark of nse state

template<typename T>
void check_in_nse(T& state, T& nse_state){
  // Input:
  //
  // state: current eos or burner state
  // nse_state: current nse state, obtained from function: get_actual_nse_state
  //
  // Purpose:
  // a boolean (state.in_nse) that tells whether it is okay to proceed with ase

  int count = 0;
  amrex::Real r = 1.0_rt;
  amrex::Real r_nse = 1.0_rt;
  
  // Check if neutron, proton, and helium are present in the network
  for (int n = 0; n < NumSpec; ++n){
    if (n == NSE_INDEX::h1_index || n == NSE_INDEX::n_index|| n == NSE_INDEX::he4_index){
      ++count;
    }
  }
  if (count < 3){
    amrex::Error("Current network does not include the proton, neutron, or helium-4, thus cannot proceed with ASE.");
  }

  // Check if n,p,a are in equilibrium
  // these two ratios are defined in the ASE paper to determine whether network is in equilibrium
  for (int n = 0; n < NumSpec; ++n){

    if (n == NSE_INDEX::h1_index || n == NSE_INDEX::n_index){
      r /= state.xn[n] * state.xn[n] * aion_inv[n] * aion_inv[n];
      r_nse /= nse_state.xn[n] * nse_state.xn[n] * aion_inv[n] * aion_inv[n];
    } 
    
    else if (n == NSE_INDEX::he4_index){
      r *= state.xn[n] * aion_inv[n];
      r_nse *= nse_state.xn[n] * aion_inv[n];
    }
  }

  std::cout << "r is " << r << std::endl;
  std::cout << "r_nse is " << r_nse << std::endl;
  
  // equilibrium condition: if pass proceed with ase if not proceed with regular eos integration
  
  if (std::abs((r - r_nse) / r_nse) < 0.5){
    state.in_nse = true;
  }
  else{
    state.in_nse = false;
  }
}

// Lets not use limiters
// template <typename T>
// void get_limiter(amrex::Real& f_full, const T& state, const amrex::Array1D<amrex::Real, 1, NumSpec+1>& ydot,
// 		 const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, const amrex::Real& t_s){

//   // Find the burning limiter
  
//   amrex::Real Y_npa = 0.0_rt;
//   amrex::Real Y_npa_dot = 0.0_rt;
//   amrex::Real Y_tilde = 0.0_rt;
//   amrex::Real Y_tilde_dot = 0.0_rt;
//   amrex::Real X_npa = 0.0_rt;
//   amrex::Real X_tilde = 0.0_rt;

//   for (int n = 0; n < NumSpec; ++n){
//     if ( n == NSE_INDEX::h1_index || n == NSE_INDEX::n_index || n == NSE_INDEX::he4_index){
//       X_npa += state.xn[n];
//       Y_npa += Y(n+1);
//       Y_npa_dot += ydot(n+1);
//     }
//     else{
//       X_tilde += state.xn[n];
//       Y_tilde += Y(n+1);
//       Y_tilde_dot += ydot(n+1);
//     }
//   }

//   // Find the burning limiters based on Equation 6 and 7 in ase paper.
//   // assume state.e is the total specific internal energy of the system.
//   std::cout << "state.e is " << state.e << std::endl;
//   std::cout << "q dot is " << ydot(NumSpec+1) << std::endl;
//   std::cout << "Y_npa_dot is " << Y_npa_dot << std::endl;
//   std::cout << "Y_tilde_dot is " << Y_tilde_dot << std::endl;
//   std::cout << "Y_npa is " << Y_npa << std::endl;
//   std::cout << "Y_tilde is " << Y_tilde << std::endl;
//   std::cout << "X_npa is " << X_npa << std::endl;
  
//   amrex::Real f_e = amrex::min(1.0_rt, state.e * f_lim / std::abs(ydot(NumSpec+1) / t_s));
//   amrex::Real f_npa = amrex::min(1.0_rt, Y_npa * f_lim / std::abs(Y_npa_dot) / t_s);
//   amrex::Real f_tilde = amrex::min(1.0_rt, Y_tilde * f_lim / std::abs(Y_tilde_dot) / t_s);

//   // f_full is the actual burning limiter we use.
//   std::cout << "f_lim is " << f_lim << std::endl;
//   std::cout << "f_e is " << f_e << std::endl;
//   std::cout << "f_npa is " << f_npa << std::endl;
//   std::cout << "f_tilde is " << f_tilde << std::endl;
  
//   f_full = amrex::min(f_e, 1.0_rt / (X_npa / f_npa + X_tilde / f_tilde));
//   std::cout << "f_full is " << f_full << std::endl;
// }


void find_max_nucs(int& max_nucs, const int current_nuc_ind){
  // Find max number of nucs possible in a reaction step given a nuclei index.
  // used to determine the size of reaction array
  
  max_nucs = 0;
  
  for (int n = NumSpec-1; n > NSE_INDEX::he4_index; --n){
    if (
	(aion[n] == aion[current_nuc_ind]-1 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]-1))
	||
	(aion[n] == aion[current_nuc_ind]+3 &&
	    (zion[n] == zion[current_nuc_ind]+1 || zion[n] == zion[current_nuc_ind]+2))
	||
	(aion[n] == aion[current_nuc_ind]+4 && zion[n] == zion[current_nuc_ind]+2)
	){
      
      ++max_nucs;
    }
    
    if (
	(aion[n] == aion[current_nuc_ind]+3 &&
	 (zion[n] == zion[current_nuc_ind]+1 || zion[n] == zion[current_nuc_ind]+2))
	 ||
	 (aion[n] == aion[current_nuc_ind]+2 &&
	  (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1
	   || zion[n] == zion[current_nuc_ind]+2))
	 ||
	 (aion[n] == aion[current_nuc_ind]-2 &&
	  (zion[n] == zion[current_nuc_ind]-2 || zion[n] == zion[current_nuc_ind]-1
	   || zion[n] == zion[current_nuc_ind]))
	){
      
      ++max_nucs;
    }
    
    if (
	(aion[n] == aion[current_nuc_ind]+2 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1
	  || zion[n] == zion[current_nuc_ind]+2))
	||
	(aion[n] == aion[current_nuc_ind]+1 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1))
	||
	(aion[n] == aion[current_nuc_ind]-3 &&
	 (zion[n] == zion[current_nuc_ind]-1 || zion[n] == zion[current_nuc_ind]-2))
	){
      
      ++max_nucs;
    }
    
    if (
	(aion[n] == aion[current_nuc_ind]+1 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1))
	||
	(aion[n] == aion[current_nuc_ind]-4 && zion[n] == zion[current_nuc_ind]-2)
	){
      
      ++max_nucs;
    }
    
  }
}

void get_rate_by_nuc(int& rate_index, amrex::Array1D<int, 1, 3> reactants_indices,
		     amrex::Array1D<int, 1, 3> products_indices){

  // fill in the corresponding rate_index by providing array of reactants and products indices
  rate_index = -1;
  
  bool found_index = false;
  
  sort_Array1D<1, 3>(reactants_indices);
  sort_Array1D<1, 3>(products_indices);
  
  for (int i = 1; i <= Rates::NumRates; ++i){
    for (int j = 1; j <= 3; ++j){
      if (NSE_INDEX::rate_indices(i, j) == reactants_indices(j) && NSE_INDEX::rate_indices(i, j+3) == products_indices(j)){
	found_index = true;
      }
      else{
	found_index = false;
      }
    }
    
    if (found_index){
      rate_index = i;
      return;
    }
    
  }

}


template <typename T>
void find_fast_reaction_cycle(const int current_nuc_ind, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
			const amrex::Array1D<amrex::Real, 1, NumSpec+1>& ydot,
			const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
			const T& state, const amrex::Real& f_full, const amrex::Real& t_s, bool& found_fast_reaction_cycle){

  // Find fast reaction cycle of a given nuclei.
  
  // Want to find max size of nucs allowed
  // template requires constant expression, but there are loops later, so it can still be used..
  
  int max_nucs;
  find_max_nucs(max_nucs, current_nuc_ind);

  // Here I use NumSpec for YHI index for now, which will leave many unused slots, but deal with it later..
  
  amrex::Array2D<int, 1, 6, 1, NumSpec+1> reactions;            // store possible reaction nuclei indices
  amrex::Array2D<int, 1, 6, 1, NumSpec+1> fast_reactions;       // store indices of fast reaction nucs

  // some variables to hold temporary data
  amrex::Array1D<int, 1, 3> reaction_scratch;
  amrex::Array1D<int, 1, 3> product_scratch;
  int rate_index;                                      // store index of the rate
  bool have_required_nucs;                             // variable used to check we have required nuc
    
  // holds forward and reverse rate of a reaction
  amrex::Real b_f = 0.0_rt;
  amrex::Real b_r = 0.0_rt;

  for (int i = 1; i <= 6; ++i){
    for (int j = 1; j <= max_nucs; ++j){
      reactions(i, j) = -1;
      fast_reactions(i, j) = -1;
    }
  }  

  // indices for reaction array to location intermediate nuclei
  
  int inter_nuc_ind1 = 1;
  int inter_nuc_ind2 = 1;
  int inter_nuc_ind3 = 1;
  int inter_nuc_ind4 = 1;

  // the beginning nuclei of a reaction cycle is always itself
  
  reactions(1, 1) = current_nuc_ind;
  
  // store all possible nuclei involved in a reaction cycle started by nuclei: current_nuc_ind
  
  for (int n = NumSpec-1; n > NSE_INDEX::he4_index; --n){
    
    if (
	(aion[n] == aion[current_nuc_ind]-1 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]-1))
	||
	(aion[n] == aion[current_nuc_ind]+3 &&
	    (zion[n] == zion[current_nuc_ind]+1 || zion[n] == zion[current_nuc_ind]+2))
	||
	(aion[n] == aion[current_nuc_ind]+4 && zion[n] == zion[current_nuc_ind]+2)
	){
      
      reactions(2, inter_nuc_ind1) = n;
      ++inter_nuc_ind1;
    }
    
    if (
	(aion[n] == aion[current_nuc_ind]+3 &&
	 (zion[n] == zion[current_nuc_ind]+1 || zion[n] == zion[current_nuc_ind]+2))
	 ||
	 (aion[n] == aion[current_nuc_ind]+2 &&
	  (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1
	   || zion[n] == zion[current_nuc_ind]+2))
	 ||
	 (aion[n] == aion[current_nuc_ind]-2 &&
	  (zion[n] == zion[current_nuc_ind]-2 || zion[n] == zion[current_nuc_ind]-1
	   || zion[n] == zion[current_nuc_ind]))
	){
      
      reactions(3, inter_nuc_ind2) = n;
      ++inter_nuc_ind2;
    }
    
    if (
	(aion[n] == aion[current_nuc_ind]+2 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1
	  || zion[n] == zion[current_nuc_ind]+2))
	||
	(aion[n] == aion[current_nuc_ind]+1 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1))
	||
	(aion[n] == aion[current_nuc_ind]-3 &&
	 (zion[n] == zion[current_nuc_ind]-1 || zion[n] == zion[current_nuc_ind]-2))
	){
      
      reactions(4, inter_nuc_ind3) = n;
      ++inter_nuc_ind3;
    }
    
    if (
	(aion[n] == aion[current_nuc_ind]+1 &&
	 (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]+1))
	||
	(aion[n] == aion[current_nuc_ind]-4 && zion[n] == zion[current_nuc_ind]-2)
	){
      
      reactions(5, inter_nuc_ind4) = n;
	++inter_nuc_ind4;
    }
  }

  reactions(5, inter_nuc_ind4) = current_nuc_ind;
  reactions(6, 1) = current_nuc_ind;

  // Checks if there are no intermediate nuclei

  for (int i = 2; i <= 4; ++i){

    have_required_nucs = false;

    for (int j = 1; j <= max_nucs; ++j){
      if (reactions(i, j) != -1){	
	have_required_nucs = true;
      }
    }
    
    if (!have_required_nucs){
      found_fast_reaction_cycle = false;

      //debug
      std::cout << "don't have required nucs" << std::endl;
      return;
    }
  }

  // print out reactions for debug for now.
  std::cout << "print out current reaction that passed init test" << std::endl;
  for (int i = 1; i <= 6; ++i){
    for (int j = 1; j <= max_nucs; ++j){
      std::cout << reactions(i, j) << " ";
    }
    std::cout << std::endl;
  }
  
  
  // The first nuc in fast reaction is always itself

  fast_reactions(1, 1) = current_nuc_ind;
  amrex::Array1D<amrex::Real, 1, 2>  Y_i;
    
  // loop over reaction stages
  for (int i = 1; i <= 5; ++i){
    // loop over nuclei in the current stage
    std::cout << "current i is " << i << std::endl;
    for (int j = 1; j <= max_nucs; ++j){
      std::cout << "current j is " << j << std::endl;
      // skip this loop if current nuc is the same as the starting nuc
      // also skip if current index is -1, which is not valid
      if (reactions(i, j) == -1 || (reactions(i, j) == current_nuc_ind && i != 1)){
	continue;
      }

      // current current stage of nuc index, (i,j) to the next stage nuc index (i+1,k)
      // get nuc index of the next stage in reactions
      
      for (int k = 1; k <= max_nucs; ++k){
	std::cout << "current k is " << k << std::endl;      
	// if reaction index is -1, skip, since its not a valid index
	if (reactions(i+1, k) == -1){
	  continue;
	}
	
	// Conditions to tell us which reaction occured.
	std::cout << "current aion is " << aion[reactions(i, j)] << std::endl;
	std::cout << "next aion is " << aion[reactions(i+1, k)] << std::endl;

	// Break following into a function,
	// take index of the reactant and product (Consistent with NumSpec)
	// need reactant_scratch, product_scratch, screened_rates, Y, state.rho, Y_i
	
	if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 4 &&
	    zion[reactions(i, j)] == zion[reactions(i+1, k)] - 2){
	  // (alpha, gamma) reaction
	  // implies reactants are he4, reactions(i, j) and product is (reactions(i+1, k)
	  // Need to get rate by providing reactant and product, get_rate_by_nuclei

	  // get indices of reaction and products:
	  // reaction indices: reactions(i, j), NSE_INDEX::he4_index
	  // product indices: reactions(i+1, k)
	  // find corresponding rate index in rate_indices

	  // fill in reaction_scratch indices

	  std::cout << "in (a,g) " << std::endl;
	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = NSE_INDEX::he4_index;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = -1;
	  product_scratch(3) = -1;
	  
	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);
	  b_f = f_full * screened_rates(rate_index) * Y(reactions(i, j) + 1) * Y(NSE_INDEX::he4_index + 1) * state.rho;
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate
	  if (rate_index == -1){
	    continue;
	  }
	  
	  b_r = f_full * screened_rates(rate_index) * Y(reactions(i+1, k) + 1);

	  Y_i(1) = Y(NSE_INDEX::he4_index + 1);
	  Y_i(2) = 0.0_rt;
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 3 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] - 2){
	  // (alpha, n) reaction

	  std::cout << " in (a,n) " << std::endl;
	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = NSE_INDEX::he4_index;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::n_index;
	  product_scratch(3) = -1;
	  
	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);
	  b_f = f_full * screened_rates(rate_index) * Y(reactions(i, j) + 1) * Y(NSE_INDEX::he4_index + 1) * state.rho;
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate
	  if (rate_index == -1){
	    continue;
	  }
	  
	  b_r = f_full * screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::n_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::he4_index + 1);
	  Y_i(2) = Y(NSE_INDEX::n_index + 1);
	}
	
	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 3 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] - 1){
	  // (alpha, p) reaction

	  std::cout << "in (a,p) " << std::endl;
	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = NSE_INDEX::he4_index;
	  reaction_scratch(3) = -1;
	  
	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::h1_index;
	  product_scratch(3) = -1;
	  
	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);
	  b_f = f_full * screened_rates(rate_index) * Y(reactions(i, j) + 1) * Y(NSE_INDEX::he4_index + 1) * state.rho;
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate
	  if (rate_index == -1){
	    continue;
	  }
	  
	  b_r = f_full * screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::h1_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::he4_index + 1);
	  Y_i(2) = Y(NSE_INDEX::h1_index + 1);
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] + 1 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)]){
	  // (gamma, n) reaction

	  std::cout << "in (g,n) " << std::endl;
	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = -1;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::n_index;
	  product_scratch(3) = -1;
	  
	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);
	  b_f = f_full * screened_rates(rate_index) * Y(reactions(i, j) + 1);
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate
	  if (rate_index == -1){
	    continue;
	  }
	  
	  b_r = f_full * screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::n_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::n_index + 1);
	  Y_i(2) = 0.0_rt;
	}
		 
	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] + 1  &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] + 1){
	  // (gamma, p) reaction

	  std::cout << "in (g,p) " << std::endl;
	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = -1;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::h1_index;
	  product_scratch(3) = -1;
	  
	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);
	  b_f = f_full * screened_rates(rate_index) * Y(reactions(i, j) + 1);
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate
	  if (rate_index == -1){
	    continue;
	  }
	  
	  b_r = f_full * screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::h1_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::h1_index + 1);
	  Y_i(2) = 0.0_rt;
	}

	else{
	  // continue if no matching reactions
	  continue;
	}

	// A condition to check whether rates are fast enough given by Eq10 and 12 in ASE paper
	std::cout << "b_f is " << b_f << std::endl;
	std::cout << "b_r is " << b_r << std::endl;
	std::cout << "before nuc is " << Y(reactions(i, j) + 1) << std::endl;
	std::cout << "after nuc is " << Y(reactions(i + 1, k) + 1) << std::endl;
	std::cout << "before nuc # is " << reactions(i, j) << std::endl;
	std::cout << "after nuc # is " << reactions(i + 1, k) << std::endl;

	std::cout << "ase_tol is " << ase_tol << std::endl;
	std::cout << "t_s is " << t_s << std::endl;


	if ((Y_i(1) / amrex::min(b_f, b_r) < ase_tol * t_s)
	    || (Y_i(2) / amrex::min(b_f, b_r) < ase_tol * t_s)
	    // && ((2.0_rt * std::abs(b_f - b_r) / (b_f + b_r)) < ase_tol)
	  ){
	
	  std::cout << "Found a fast reaction" << std::endl;
	  // store nuclei index to fast_reaction if pass the condition
	  std::cout << "Now add nuc " << reactions(i+1, k) << " to fast reaction " << std::endl;
	  fast_reactions(i + 1, k) = reactions(i + 1, k);
	}
     		      
      }
    }
    // if (fast_reactions(i + 1, p) have any non -1, means found a fast reac nuc
    
    have_required_nucs = false;
    
    for (int p = 1; p <= max_nucs; ++p){

      // replace reaction nucs with fast_reaciton nucs of the next reaction stage

	reactions(i + 1, p) = fast_reactions(i + 1, p);

	// if we see current_nuc_ind in the next stage, means we have finished the cycle

	if (fast_reactions(i + 1, p) == current_nuc_ind){
	  found_fast_reaction_cycle = true;
	  return;
	}

	if (fast_reactions(i + 1, p) != -1){
	  have_required_nucs = true;
	}
    }

    // check if we even have nucs in the next stage that are fast enough
    // if not, return
    
    if (!have_required_nucs){
      found_fast_reaction_cycle = false;
      std::cout << "did not find fast reaction in the next stage, breaking in current nuc " << std::endl;
      return;
    }
  }
    
}


template <typename T>
void actual_check_in_nse(T& state, T& nse_state){
  
  // Runtime parameters:
  // f: characteritic scale for the burning limiter, ~0.05
  // cell_dx: size of simulation mesh, ~1.0e6
  // 
  // assume we have f, cell_dx in _parameters
  
  // Update whether state is in the ballpark of NSE
  check_in_nse(state, nse_state);
  if (!state.in_nse){
    std::cout << "state currently not in NSE, switch to eos!" << std::endl;
    return;
  }

  // set molar fractions
  // initialize ydot to store Ydot and energy generation rate.
  
  amrex::Array1D<amrex::Real, 1, NumSpec> Y;
  amrex::Array1D<amrex::Real, 1, NumSpec + 1> ydot;
  
  for (int n = 1; n <= NumSpec; ++n){
    Y(n) = state.xn[n-1] * aion_inv[n-1];
    ydot(n) = 0.0_rt;
  }


  // Find ydot and energy generation rate
  
#ifdef NEW_NETWORK_IMPLEMENTATION
  // rhs_t rhs_data(rate), rate is a string (enum to be indices);
  // gives specie names, which can be used to determine reaction and find the corresponding rate
  amrex::Error("Not implemented yet");
  //constexpr rhs_t data = RHS::rhs_data(rate)
  // Real rhs_term<int species, int rate> (state, Real fr, Real rr), given species, rate, forward and reverse.
  
  //RHS::rhs(state, ydot);
  //actual_rhs(state, ydot);
#else
    
  // Get rates, ydot,and energy generation rate, copy routine from  rhs_nuc(state, ydot)
  
  // Need individual rates in rate_eval
  rate_t rate_eval;
  constexpr int do_T_derivatives = 0;
  evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

  eos_t eos_state;                             // need eos_state for speed of sound
  amrex::Real t_s;                             // a parameter to characterize whether a rate is fast enough

  if constexpr (std::is_same<T, eos_t>::value){
    t_s = cell_dx / state.cs;           
  }
  else{
    burn_to_eos<eos_t>(state, eos_state);
    eos(eos_input_rt, eos_state);
    t_s = cell_dx / eos_state.cs;
  }
  
  // Find ydot
  rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

  // Find energy generation rate
  amrex::Real enuc;
  ener_gener_rate(ydot, enuc);

#ifdef NEUTRINOS
  amrex::Real sneut, dsneutdt, dsneutdd, snuda, snudz;
  sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);
#else
  amrex::Real sneut = 0.0_rt, dsneutdt =0.0_rt , dsneutdd = 0.0_rt, snuda = 0.0_rt, snudz = 0.0_rt;
#endif
  // fill in energy generation rate to ydot
  ydot(NumSpec + 1) = enuc - sneut;
#endif

  // Limiter is not desired for now, so lets set them to 1 for now.
  amrex::Real f_full = 1.0_rt;
  // get_limiter(f_full, state, ydot, Y, t_s);

  
  // Now we look through the network and see if there are fast reaction cycles
  // Need to separate forward and reverse rate and detemine each step is fast enough.
  // use vectors for now
  
  bool found_fast_reaction_cycle = false;
  
  // Do a reverse for loop to start from heaviest nuclei
  
  for (int n = NumSpec-1; n >= 0; --n){
    if (found_fast_reaction_cycle){
      break;
    }

    find_fast_reaction_cycle(n, Y, ydot, rate_eval.screened_rates, state, f_full, t_s, found_fast_reaction_cycle);
  }

  // Change state.in_nse accordingly
  // Do nse grouping if found fast reaction cycle.

  amrex::Array1D<int, 1, NumSpec> group_ind;        // holds group index for each nuclei
  
  if (!found_fast_reaction_cycle){
    std::cout << "Unable to find fast reaction cycle, so cannot proceed with ase" << std::endl;
    state.in_nse = false;
  }
  else{
    // Do nse grouping 
    // state.in_nse = true;
    std::cout << "Found fast reaction" << std::endl;

    nse_grouping(group_ind, state, Y, rate_eval.screened_rates, t_s);
    
    // Check if we result in a single group:
    if (in_single_group<1, NumSpec>(group_ind)){
      state.in_nse = true;
      std::cout << "In a single group" << std::endl;
    }
    else{
      state.in_nse = false;
      //debug
      for (int z = 1; z <= NumSpec; ++z){
	std::cout << group_ind(z) << std::endl;
      }
      std::cout << "Not in a single group" << std::endl;
    }
  }
}


// After all checks are successful, lets do nse grouping.

int get_root_index(const int nuc_ind,
		     const amrex::Array1D<int, 1, NumSpec>& group_ind){

  // This function returns the root index of the nuclei
  // by providing the nuclei index [0, NumSpec-1], and group indices, group_ind

  int root_index;
  int scratch_ind = nuc_ind;
  
  while(true){
    root_index = group_ind(scratch_ind + 1);

    if (root_index != scratch_ind + 1){
      scratch_ind = root_index - 1;
    }

    else{
      return root_index;
    }
  }
  
}


void nse_union(const int nuc_ind_a, const int nuc_ind_b, amrex::Array1D<int, 1, NumSpec>& group_ind){

  // This function joins the two group of the two nuc indices:nuc_ind_a and nuc_ind_b
  // The smaller group is joined to the larger group.

  int root_index_a = get_root_index(nuc_ind_a, group_ind);;
  int root_index_b = get_root_index(nuc_ind_b, group_ind);;

  if (root_index_a == root_index_b){
    return;
  }
  
  // find size of the two groups containing a and b
  
  int group_a_size = 0;
  int group_b_size = 0;

  for (int n = 1; n <= NumSpec; ++n){
    if (n == root_index_a){
      ++group_a_size;
    }
    if (n == root_index_b){
      ++group_b_size;
    }
  }

  // If there is only 1 isotope that is not in LIG, merge it into LIG
  if (group_a_size >= group_b_size){
    group_ind(root_index_b) = group_ind(root_index_a);
  }
  
  else{
    group_ind(root_index_a) = group_ind(root_index_b);
  }
}


void fill_pair_rate_index(int& pair_rate_index, const int& rate_index){

  // This funciton finds the pair rate index given a rate.
  // given forward rate index, return reverse rate index, vice versa

  pair_rate_index = -1;

  for (int n = 1; n <= Rates::NumRates; ++n){
    for (int i = 1; i <= 3; ++i){
      if ((NSE_INDEX::rate_indices(n, i) == NSE_INDEX::rate_indices(rate_index, i+3))
	  && (NSE_INDEX::rate_indices(n, i+3) == NSE_INDEX::rate_indices(rate_index, i))
	  ){
	pair_rate_index = n;
      }
    }
    
  }
}


// void add_eliminated_rate(int& index, const int current_rate,
// 			 amrex::Array1D<int, 1, Rates::NumRates>& eliminated_rates,
// 			 const amrex::Array1D<int, 1, NumSpec>& group_ind){
  

//   // this function adds current rate to eliminated rate based on conditions

//   // check if current_rate_index is already in eliminated_rates

//   for (int n = 1; n <= Rates::NumRates; ++n){
//     if (index == eliminated_rates(n)){
//       return;
//     }
//   }

//   // if there are no reverse rate, then add that rate to eliminated_rates
//   int pair_rate_index;
//   int num_nonLIG = 0;

//   // find the reverse rate index
  
//   fill_pair_rate_index(pair_rate_index, current_rate);

//   // Check if there are > 2 isotopes or 0 isotopes in LIG
    
//   for (int k = 1; k <= 6; ++k){

//     if (NSE_INDEX::rate_indices(n, k) == -1){
//       continue;
//     }
      
//     root_index = get_root_index(NSE_INDEX::rate_indices(current_rate, k), group_ind);

//     if (root_index != 1){
//       ++num_nonLIG;
//     }
//   }

//   // 3 conditions:
//   //
//   // 1) add rates to eliminated rates when reaction has more than 3 reactant or product
//   // 2) there are 3 reactants or products involved in the network
//   // 3) number of nonLIG isotopes are more than 2 or none.
  
//   if ( (pair_rate_index == -1)
//        || (NSE_INDEX::rate_indices(current_rate, 1) != -1)
//        || (NSE_INDEX::rate_indices(current_rate, 4) != -1)
//        || (num_nonLIG > 2)
//        || (num_nonLIG == 0)
//       ){
//     eliminated_rates(index) = current_rate;
//     ++index;
//   }
  
// }



template<typename T>
void is_fastest_rate(amrex::Array1D<int, 1, 2>& merge_indices,
		     amrex::Real& rate_t, const int current_rate,
		     const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, const T& state,
		     const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
		     const amrex::Array1D<int, 1, NumSpec>& group_ind,
		     const amrex::Real& t_s){
  
  // This function determines whether current_rate has a faster timescale than
  // the provided time scale, rate_t.
  // It also skips current_rate if certain condition is met.

  int pair_rate_index;
  int root_index;
  int num_nonLIG = 0;
  int nonLIG_root = -1;
  
  // Find the Y_group of the nuclei involved in the net

  amrex::Array1D<amrex::Real, 1, 6> Y_group;

  // Here Y_nonNPA represent 2 group isotope that are not n,p,a

  amrex::Array1D<int, 1, 2> nonNPA_ind = {-1, -1};
  int m;
  
  // find the reverse rate index
  
  fill_pair_rate_index(pair_rate_index, current_rate);

  // Check if there are > 2 isotopes or 0 isotopes in LIG
    
  for (int k = 1; k <= 6; ++k){

    if (NSE_INDEX::rate_indices(current_rate, k) == -1){
      continue;
    }
      
    root_index = get_root_index(NSE_INDEX::rate_indices(current_rate, k), group_ind);

    // Determine number of nonLIG isotopes
    // also check whether nonLIG isotopes are already merged
    
    if (root_index != 1){
      ++num_nonLIG;

      if (root_index != nonLIG_root){
	nonLIG_root = root_index;
      }

      // return if nonLIG_root index is repeated, meaning these isotope already merged
      else{
	return;
      }
    }

    // Find the group molar fraction of the current isotope
    
    for (int n = 0; n < NumSpec; ++n){
	  
      // if nuclei have the same root_index, then they are the same group
      if (get_root_index(n, group_ind) == root_index){	    
	Y_group(k) += Y(n + 1);
      }
    }

    // Find indices that is non n,p,a
    
    if (NSE_INDEX::rate_indices(current_rate, k) != NSE_INDEX::h1_index
	|| NSE_INDEX::rate_indices(current_rate, k) != NSE_INDEX::n_index
	|| NSE_INDEX::rate_indices(current_rate, k) != NSE_INDEX::he4_index
	){

      // return if there are more than 2 nonNPA nuclei involved in the network
      
      if (m > 2){
	return;
      }
      
      nonNPA_ind(m) = k;
      ++m;
    }
    
  }

  // 3 other conditions to skip current rate:
  //
  // 1) skip if there is no reverse rate involved.
  // 2) add rates to eliminated rates when reaction has more than 3 reactant or product
  // 3) number of nonLIG isotopes are more than 2 or none.
  
  if ( (pair_rate_index == -1)
       || (NSE_INDEX::rate_indices(current_rate, 1) != -1)
       || (NSE_INDEX::rate_indices(current_rate, 4) != -1)
       || (num_nonLIG > 2)
       || (num_nonLIG == 0)
      ){
    return;
  }
   
  // Find the rates
  
  amrex::Real b_f;
  amrex::Real b_r;
  
  b_f = screened_rates(current_rate) * Y_group(3);  
  b_r = screened_rates(pair_rate_index) * Y_group(6);

  if (NSE_INDEX::rate_indices(current_rate, 2) != -1){
    b_f *= Y_group(2) * state.rho;
  }

  if (NSE_INDEX::rate_indices(current_rate, 5) != -1){
    b_r *= Y_group(5) * state.rho;
  }

  // Find the timescale of the rate
  amrex::Real scratch_t = amrex::min(Y(nonNPA_ind(1)), Y(nonNPA_ind(2))) / amrex::min(b_f, b_r);

  // return if current rate time scale is larger than previous time scale
  if (rate_t < scratch_t){
    return;
  }
  
  if (scratch_t < ase_tol * t_s
      && 2.0_rt * std::abs(b_f - b_r) / (b_f + b_r) < ase_tol
      ){

    rate_t = scratch_t;    
    merge_indices(1) = NSE_INDEX::rate_indices(current_rate, nonNPA_ind(1));
    merge_indices(2) = NSE_INDEX::rate_indices(current_rate, nonNPA_ind(2));
  }
  
}



template<typename T>
void nse_grouping(amrex::Array1D<int, 1, NumSpec>& group_ind, const T& state,
		  const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
		  const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
		  const amrex::Real& t_s){
  
    
  // Make sure check_in_nse() is called before nse_grouping

  if (!state.in_nse){
    amrex::Error("State is not in the ballpark of NSE, call check_in_nse() before this");
  }
  
  // fill in initial group_ind, group_ind go from 1 to NumSpec

  for (int i = 1; i <= NumSpec; ++i){
    group_ind(i) = i;
    
    // let n,p,a form the same group (LIG) initially, let 1 be index of LIG
    if (i == NSE_INDEX::h1_index + 1 || i == NSE_INDEX::n_index + 1
	|| i == NSE_INDEX::he4_index + 1){
      group_ind(i) = 1;
    }
  }

  amrex::Array1D<int, 1, 2> merge_indices;
  amrex::Real rate_t = std::numeric_limits<amrex::Real>::max();

  // amrex::Array1D<int, 1, Rates::NumRates> eliminated_rates;
  // int eliminated_rate_index = 1;
  
  // // initial process to eliminated rates
  // for (int n = 1; n <= Rates::NumRates; ++n){

  //   add_eliminated_rate(eliminated_rate_index, n, eliminated_rates, group_ind);
  // }
  
  
  // Go through each reaction and perform grouping
  while(true){
    
    merge_indices(1) = -1;
    merge_indices(2) = -1;

    // Find the fastest time scale of the avaliable reaction
    for (int n = 1; n <= Rates::NumRates; ++n){

      // Check if current rate is the faster rate than previous rate
      is_fastest_rate(merge_indices, rate_t, n, Y, state, screened_rates, group_ind, t_s);
    }

    // if there are no reaction satisfy conditions, break
    
    if (merge_indices(1) == -1 || merge_indices(2) == -1){
      break;
    }

    // union the isotopes of the rate that have the fastest time scale
    nse_union(merge_indices(1), merge_indices(2), group_ind);
  }
  
}

#endif
