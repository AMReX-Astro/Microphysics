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


void cap_first_letter(std::string& name){
  // capitalize first letter of nuc except for neutron and proton, n and p.
  
  if (name == "p" || name == "n"){
    return;
  }
  
  name[0] = toupper(name[0]);
}


void find_rate_pair(const std::string& rate_name, std::string& reverse_rate,
		    std::vector<std::string>& reactants, std::vector<std::string>& products){
  // find the reverse rate_name, reactants and products in a given rate, Rates::rate_names
  // store them vector of strings.
  //
  // currently not used
  //
  
  reverse_rate = "";
  reactants = {};
  products = {};
  
  auto pos = rate_name.find("__");
  std::string reactant_names = rate_name.substr(0, pos);
  std::string product_names = rate_name.substr(pos+2);
  
  pos = product_names.find("__");
  std::string suffix = product_names.substr(pos);
  product_names.erase(pos);
  
  std::string reactant;
  std::string product;
  
  // store names of reactants and products

  while ((pos = product_names.find('_')) != std::string::npos){
    product = product_names.substr(0, pos);

    reverse_rate += product;
    reverse_rate += "_";
    
    cap_first_letter(product);
    products.push_back(product);
    product_names.erase(0, pos+1);
  }

  reverse_rate += product_names;
  reverse_rate += "__";
  
  cap_first_letter(product_names);
  products.push_back(product_names); 
  
  while ((pos = reactant_names.find('_')) != std::string::npos){
    reactant = reactant_names.substr(0, pos);

    reverse_rate += reactant;
    reverse_rate += "_";
    
    cap_first_letter(reactant);
    reactants.push_back(reactant);
    reactant_names.erase(0, pos+1);
  }
  reverse_rate += reactant_names;
  reverse_rate += suffix;
  
  cap_first_letter(reactant_names);
  reactants.push_back(reactant_names);
}


void find_rate_indices(amrex::Array1D<int, 1, 3>& reaction_indices,
		       amrex::Array1D<int, 1, 3>& product_indices, const std::string& rate_name){

  // fill the corresponding reactant and product indices of a given rate
  
  for (int n = 1; n <= 3; ++n){
    reaction_indices(n) = -1;
    product_indices(n) = -1;
  }

  int i;

  // Process rate_name to find the reactants and products
  
  auto pos = rate_name.find("__");
  std::string reactant_names = rate_name.substr(0, pos);
  std::string product_names = rate_name.substr(pos+2);
  
  pos = product_names.find("__");
  std::string suffix = product_names.substr(pos);
  product_names.erase(pos);

  std::string _temp;

  // Find the index of reactants and products in short_spec_names
  i = 1;
  while ((pos = reactant_names.find('_')) != std::string::npos){

    _temp = reactant_names.substr(0, pos);

    // aprox network need capitalize
    // but we only work with pynucastro networks..
  
    // cap_first_letter(reactant);    

    for (int n = 0; n < NumSpec; ++n){
      if (short_spec_names_cxx[n] == _temp){
	reaction_indices(i) = n;
	++i;
      }
    }
        
    reactant_names.erase(0, pos+1);
  }
   
  // cap_first_letter(reactant_names);
  for (int n = 0; n < NumSpec; ++n){
    if (short_spec_names_cxx[n] == reactant_names){
      reaction_indices(i) = n;
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
 
  sort_Array1D<1, 3>(reaction_indices);
  sort_Array1D<1, 3>(product_indices);
  
}

template <int l, int m>
void sort_Arrary1D(amrex::Array1D<int, l, m>& array_1d){

  // a simple bubble sort, for sorting Array1D
  // assumes index start from int l, and end at int m

  bool sorted = false;
  int scratch;
  
  for (int i = l; i <= m-1; ++i){
    for (int j = l; j <= m-i-1; ++j){
      if (array_1d(j) > array_1d(j+1)){
	scratch = array_1d(j);
	array_1d(j) = array_1d(j+1);
	array_1d(j+1) = scratch;
      }
    }
  }
}


void init_actual_nse(){
  
  // Initialization for nse, check whether network is valid and stores indices
  
  for (int n = 0; n < NumSpec; ++n){

    // store index of photoionization proton by looking for Z=1 but
    // exclude h1 which is always at index=0

    if (zion[n] == 1 && aion[n] == 1){
      if (n == 0){
	NSE_INDEX::h1_index = n;
      }
      else {
	NSE_INDEX::p_index = n;
      }
    }

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

  amrex::Array1D<int, 1, 3> reaction_indices;
  amrex::Array1D<int, 1, 3> product_indices;
  
  for (int n = 1; n <= Rates::NumRates; ++n){

    // get current name of the rate
    current_rate_name = Rates::rate_name[n];

    // find the reaction_indices and product indices
    find_rate_indices(reaction_indices, product_indices, current_rate_name);

    // fill in rate_indices
    for (int i = 1; i <= 3; ++i){
      rate_indices(n, i) = reaction_indices(i);
      rate_indices(n, i+3) = product_indices(i);
    }
    
    rate_indices(n, 7) = n;    
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


// Check whether to proceed with Adaptive statistical equilibrium
template<typename T>
void ase_check(T& state, T& nse_state){
  // Input:
  //
  // state: current eos or burner state
  // nse_state: current nse state, obtained from function: get_actual_nse_state
  //
  // Purpose:
  // a boolean (state.in_nse) that tells whether it is okay to proceed with ase

  int count{0};
  amrex::Real r = 1.0_rt;
  amrex::Real r_nse = 1.0_rt;
  
  // Check if neutron, proton, and helium are present in the network
  for (int n = 0; n < NumSpec; ++n){
    if (NSE_INDEX::h1_index == n || NSE_INDEX::n_index == n || NSE_INDEX::he4_index == n){
      ++count;
    }
}
  if (count < 3){
    amrex::Error("Current network does not include the proton, neutron, or helium-4, thus cannot proceed with ASE.");
  }

  // Check if n,p,a are in equilibrium
  // these two ratios are defined in the ASE paper to determine whether network is in equilibrium
  for (int n = 0; n < NumSpec; ++n){

    if (NSE_INDEX::h1_index == n || NSE_INDEX::n_index = n){
      r /= std::sqrt(state.xn[n] * aion_inv[n]);
      r_nse /= std::sqrt(nse_state.xn[n] * aion_inv[n]);
    } 
    
    else if (NSE_INDEX::h4_index == n){
      r *= state.xn[n] * aion_inv[n];
      r_nse *= nse_state.xn[n] * aion_inv[n];
    }
  }

  // equilibrium condition: if pass proceed with ase if not proceed with regular eos integration
  if (std::abs((r - r_nse) / r) < 0.5){
    state.in_nse = true;
  }
  else{
    state.in_nse = false;
  }
}


template <typename T>
void get_limiter(amrex::Real& f_full, const T& state, const amrex::Array1D<amrex::Real, 1, NumSpec+1>& ydot,
		 const amrex::Array1D<amrex::Real, 1, NumSpec>& Y){

  // Find the burning limiter
  
  amrex::Real Y_npa = 0.0_rt;
  amrex::Real Y_npa_dot = 0.0_rt;
  amrex::Real Y_tilde = 0.0_rt;
  amrex::Real Y_tilde_dot = 0.0_rt;
  amrex::Real X_npa = 0.0_rt;
  amrex::Real X_tilde = 0.0_rt;

  for (int n = 0; n < NumSpec; ++n){
    if (NSE_INDEX::p_index == n){
      continue;
    }

    if (NSE_INDEX::h1_index == n || NSE_INDEX::n_index == n || NSE_INDEX::he4_index == n){
      X_npa += state.xn[n];
      Y_npa += Y(n+1);
      Y_npa_dot += ydot(n+1)
    }
    else{
      X_tilde += state.xn[n];
      Y_tilde += Y(n+1);
      Y_tilde_dot += ydot(n+1)
    }
  }

  // Find the burning limiters based on Equation 6 and 7 in ase paper.
  // assume state.e is the total specific internal energy of the system.

  amrex::Real f_e = amrex::min(1.0_rt, state.e * f_lim / std::abs(enuc) / t_s);
  amrex::Real f_npa = amrex::min(1.0_rt, Y_npa * f_lim / std::abs(Y_npa_dot) / t_s);
  amrex::Real f_tilde = amrex::min(1.0_rt, Y_tilde * f_lim / std::abs(Y_tilde_dot) / t_s);

  // f_full is the actual burning limiter we use.
  
  f_full = amrex::min(f_e, 1.0_rt / (X_npa / f_npa + X_tilde / f_tilde));
}


void find_max_nucs(int& max_nucs, const int current_nuc_ind){
  // Find max number of nucs possible in a reaction step given a nuclei index.
  // used to determine the size of reaction array
  
  max_nucs = 0;
  
  for (int n = NumSpec-1, n > NSE_INDEX::he4_index; --n){

    if (aion[n] == aion[current_nuc_ind]-1 &&
	(zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]-1)){

      ++max_nucs;
    }

    else if (aion[n] == aion[current_nuc_ind]-2 && (zion[n] == zion[current_nuc_ind]
       || zion[n] == zion[current_nuc_ind]-1 || zion[n] == zion[current_nuc_ind]-2)){

      ++max_nucs;
    }
    
    else if (aion[n] == aion[current_nuc_ind]-3 &&
	(zion[n] == zion[current_nuc_ind]-1 || zion[n] == zion[current_nuc_ind]-2)){
      
      ++max_nucs;
    }
    
    else if (aion[n] == aion[current_nuc_ind]-4 && zion[n] == zion[current_nuc_ind]-2){

      ++max_nucs;
    }
  }
}

void find_fast_reaction(const int current_nuc_ind){
  // Find fast reaction of a given nuclei.

  // find maximum nucs possible to declare size of reactions
  
  int max_nucs;
  find_max_nucs(max_nucs, current_nuc_ind);
  
  amrex::Array2D<int, 1, 6, 1, max_nucs> reactions;            // store possible reaction nuclei indices

  for (int i = 1; i <= 6; ++i){
    for (int j = 1; j <= max_nucs; ++j){
      reactions(i, j) = -1;
    }
  }  
  
  int inter_nuc_ind1 = 1;
  int inter_nuc_ind2 = 1;
  int inter_nuc_ind3 = 1;
  int inter_nuc_ind4 = 1;

  // the beginning nuclei of a reaction cycle is always itself
  
  reactions(1, 1) = current_nuc_ind;
  
  // store all possible nuclei involved in a reaction cycle started by nuclei: current_nuc_ind
  
  for (int n = NumSpec-1, n > NSE_INDEX::he4_index; --n){

      if (aion[n] == aion[current_nuc_ind]-1 &&
	  (zion[n] == zion[current_nuc_ind] || zion[n] == zion[current_nuc_ind]-1)){

	reactions(2, inter_nuc_ind1) = n;
	++inter_nuc_ind1;
      }

      if (aion[n] == aion[current_nuc_ind]-2 && (zion[n] == zion[current_nuc_ind]
	     || zion[n] == zion[current_nuc_ind]-1 || zion[n] == zion[current_nuc_ind]-2)){

	reactions(3, inter_nuc_ind2) = n;
	++inter_nuc_ind2;
      }

      if (aion[n] == aion[current_nuc_ind]-3 &&
	  (zion[n] == zion[current_nuc_ind]-1 || zion[n] == zion[current_nuc_ind]-2)){

	reactions(4, inter_nuc_ind3) = n;
	++inter_nuc_ind3;
      }

      if (aion[n] == aion[current_nuc_ind]-4 && zion[n] == zion[current_nuc_ind]-2){

	reactions(5, inter_nuc_ind4) = n;
	++inter_nuc_ind4;
      }
    }

  reactions(5, inter_nuc_ind4) = current_nuc_ind;
  reactions(6, 1) = current_nuc_ind;

  // Checks if there are no intermediate nuclei

  bool have_intermediate_nucs; 
  for (int i = 2; i <= 4; ++i){

    have_intermediate_nucs = false;

    for (int j = 1; j <= max_nucs; ++j){
      if (reactions(i, j) != -1){
	have_intermediate_nucs = true;
      }
    }
    
    if (!have_intermediate_nucs){
      return;
    }
  }
  

  // determine where the reactions are fast reactions

  amrex::Array2D<int, 1, 6, 1, max_nucs> fast_reactions;
  amrex::Array1D<int, 1, 3> reaction_scratch;
  amrex::Array1D<int, 1, 3> product_scratch;

  // initialize fast_reactions indices
  
  for (int i = 1; i <= 6; ++i){
    for (int j = 1; j <= max_nucs; ++j){
      fast_reactions(i, j) = -1;
    }
  }  

  fast_reactions(1, 1) = current_nuc_ind;

  // loop over reaction stages
  for (int i = 1; i <= 5; ++i){
    // loop over nuclei in the current stage
    for (int j = 1; j <= max_nucs; ++j){

      // skip this loop if current nuc is the same as the starting nuc
      if (reactions(i, j) == current_nuc_ind && i != 1){
	continue;
      }

      // current current stage of nuc index, (i,j) to the next stage nuc index (i+1,k)
      // get nuc index of the next stage in reactions
      
      for (int k = 1; k <= max_nucs; ++k){

	// Conditions to tell us which reaction occured.
	
	if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 4){
	  // (alpha, gamma) reaction
	  // implies reactants are he4, reactions(i, j) and product is (reactions(i+1, k)
	  // Need to get rate by providing reactant and product, get_rate_by_nuclei

	  // get indices of reaction and products:
	  // reaction indices: reactions(i, j), NSE_INDEX::he4_index
	  // product indices: reactions(i+1, k)
	  // find corresponding rate index in rate_indices
	  
	  
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 3 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] - 2){
	  // (alpha, n) reaction
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 3 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] - 1){
	  // (alpha, p) reaction
	}

	else if (zion[reactions(i, j)] == zion[reactions(i+1, k)]){
	  // (n, gamma) reaction
	}

	else if (zion[reactions(i, j)] == zion[reactions(i+1, k)] + 1){
	  // (p, gamma) reaction
	}

	else{
	  // continue if no matching reactions
	  
	  continue;
	}

	// Find the forward and reverse rate of the given reaction

	// A condition to check whether rates are fast enough given by Eq10 and 12 in ASE paper

	// if there are fast reaction nuclei found in the next stage, replace nuclei in reactions with those fast nuclei
	
      }
      
    }
  }
    
}





template <typename T>
void ase (T& state, T& nse_state){

  // Runtime parameters:
  // f: characteritic scale for the burning limiter, ~0.05
  // cell_dx: size of simulation mesh, ~1.0e6
  // 
  // assume we have f, cell_dx in _parameters
  
  // Update whether state is about in NSE in order to do ASE.
  ase_check(state, nse_state);
  if (!state.in_nse){
    amrex::Error("state currently not in NSE, switch to eos!");
  }

  // set molar fractions
  // initialize ydot to store Ydot and energy generation rate.
  
  amrex::Array1D<amrex::Real, 1, NumSpec> Y;
  amrex::Array1D<amrex::Real, 1, NumSpec + 1> ydot;
  for (int n = 1; n <= NumSpec; ++n){
    if (NSE_INDEX::p_index == n-1){
      continue;
    }
    Y(n) = state.xn[n-1] * inv_aion[n-1];
    ydot(n) = 0.0_rt;
  }

  // Need speed of sound and simulation mesh size, make dx: global variable ~ 1.0e6
  
  amrex::Real t_s = cell_dx / state.cs;

  // Find ydot and energy generation rate
  // Currently also considering neutrino thermal loss if neutrino included
  
#ifdef NEW_NETWORK_IMPLEMENTATION
  // rhs_t rhs_data(rate), rate is a string (enum to be indices);
  // gives specie names, which can be used to determine reaction and find the corresponding rate
  amrex::error("Not implemented yet");
  //constexpr rhs_t data = RHS::rhs_data(rate)
  // Real rhs_term<int species, int rate> (state, Real fr, Real rr), given species, rate, forward and reverse.
  
  //RHS::rhs(state, ydot);
  //actual_rhs(state, ydot);
#else
    
  // Get rates, ydot,and energy generation rate, copy routine from  rhs_nuc(state, ydot)
  
  // Need individual rates in rate_eval
  rate_t rate_eval;
  constexpr do_T_derivatives = 0;
  evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);
  
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
  
  ydot(NumSpec + 1) = enuc - sneut;
  //  actual_rhs(state, ydot);
#endif

  // Get burning limiter
  
  amrex::Real f_full;
  get_limiter(f_full, state, ydot, Y);
  
  // Now we look through the network and see if there are fast reaction cycles
  // Need to separate forward and reverse rate and detemine each step is fast enough.
  // use vectors for now
  
  bool found_fast_reaction = false;
  
  // Do a reverse for loop to start from heaviest nuclei
  
  for (int n = NumSpec-1; n => 0; --n){
    if (found_fast_reaction){
      break;
    }
        
    find_fast_reactions(n);
    
  }

}


#endif




  
  // Strip all rate names, to those include p, n or he4

  
  // A procedure to find the pair rates of a given rate_name
  amrex::Real b_f;
  amrex::Real b_r;

  // Assume we need all intermediate nuclei
  // only need b_f and b_r so can get below into a same routine and find b_f and b_r

  // Should keep track of what rates are used, put used rates into a list.
  std::string f_rate_name;
  std::string r_rate_name;
  std::vector<std::string> reactants;
  std::vector<std::string> products;
  
  for (int i = 1; i <= Rates::NumRates; ++i){
    f_rate_name = Rates::rate_names[i];
    find_rate_pair(f_rate_name, r_rate_name, reactants, products);
    
    b_f = rate_eval.screened_rates(i);

    // Find corresponding reverse rate
    // uses short_spec_name_cxx for now...
    for (int n = 1; n <= Rates::NumRates; ++n){
      if (rate_names[n] == r_rate_name){
	b_r = rate_eval.screened_rates(n);
      }
    }

    // Need molar fractions of reactants
    for (int n = 0; n < NumSpec; ++n){
      for (auto nuc : reactants){
	if (short_spec_name_cxx[n] == nuc){
	  b_f *= Y(n+1);
	}
      }

      for (auto nuc : products){
	if (short_spec_names_cxx[n] == nuc){
	  b_r *= Y(n+1);
	}
      }
    }

    b_f *= f_full;
    b_r *= f_full;
  }
