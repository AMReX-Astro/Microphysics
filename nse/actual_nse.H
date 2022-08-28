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


namespace NSE_INDEX
{
  extern bool initialized;
  
  // store indices:
  extern AMREX_GPU_MANAGED int p_index;
  
  // these three indices are helpful for implementing adaptive statistical equilibrium.  
  extern AMREX_GPU_MANAGED int h1_index;
  extern AMREX_GPU_MANAGED int n_index;
  extern AMREX_GPU_MANAGED int he4_index;
}


// Initialization for nse, check whether network is valid and stores indices
void init_actual_nse(){
  
  for (n = 0; n < NumSpec; ++n){

    // store index of photoionization proton by looking for Z=1 but exclude h1 which is always at index=0
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
      NSE_INDEX::he_index = n;
    }
  }

  
  // Check if network results in singular jacobian first, require at least one nuclei that nuclei.Z != nuclei.N
  // Some examples include aprox13 and iso7
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

  NSE_INDEX::initialized = true;
}
  

// struct to store output from constraint equations and jacobian for newton-raphson
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
    gamma = std::pow(zion[n], 5.0_rt/3.0_rt) * C::q_e * C::q_e * std::cbrt(4.0_rt * M_PI * n_e / 3.0_rt) / C::k_B / state.T;

    // chemical potential for coulomb correction
    u_c = C::k_B * state.T / C::Legacy::MeV2erg * (A1 * (std::sqrt(gamma * (A2 + gamma)) - A2 * std::log(std::sqrt(gamma / A2) + std::sqrt(1.0_rt + gamma / A2))) + 2.0_rt * A3 * (std::sqrt(gamma) - std::atan(std::sqrt(gamma))));

    // find nse mass frac
    nse_state.xn[n] = network::mion(n+1) * partition_function / state.rho * std::pow(2.0 * M_PI * network::mion(n+1) * C::k_B * state.T / std::pow(C::hplanck, 2.0_rt), 3.0_rt/2.0_rt) *  std::exp((zion[n] * state.mu_p + (aion[n] - zion[n]) * state.mu_n - u_c + network::bion(n+1)) / C::k_B / state.T * C::Legacy::MeV2erg);
  }

  // assign temperature and density
  nse_state.T = state.T;
  nse_state.rho = state.rho;

  return nse_state;
}

// get the constraint equations and its jacobian for newton raphson
template <typename T>
Newton_inputs nse_constraint(const T& state){
  // This functions finds the constraint equations and jacobian used for calculating nse.

  // it is used to store constraint equation and jacobian
  Newton_inputs nse_inputs;

  // calculate the nse state based on initial conditions of mu_p and mu_n
  // which are chemical potential of proton and neutron
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
  nse_inputs.eqs(1) =  nse_state.y_e - state.y_e;

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

// A newton-raphson solver for finding nse state used for calibrating chemical potential of proton and neutron
template<typename T>
void nse_nr_solver(T& state, amrex::Real eps=1.0e-3_rt) {

  bool converged = false;                                     // whether nse solver converged or not

  Newton_inputs f = nse_constraint(state);                    // get constraint eqs and jacobian

  amrex::Real det;                                            // store determinant for finding inverse jac
  decltype(f.jac) inverse_jac;                                // store inverse jacobian
  amrex::Real d_mu_p = 0.0_rt;                                // difference in chemical potential of proton
  amrex::Real d_mu_n = 0.0_rt;                                // difference in chemical potential of neutron

  // begin newton-raphson
  for (int i = 0; i < max_nse_iters; ++i){

    // check if current state fulfills constraint equation
    if (std::abs(f.eqs(0)) < eps && std::abs(f.eqs(1)) < eps){
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

  
  // Check whether input ye is actually valid
  amrex::Real ye_lo = 1.0_rt;
  amrex::Real ye_hi = 0.0_rt;
  
  for (int n = 0; n < NumSpec; ++n){
    ye_lo = amrex::min(zion[n] * aion_inv[n], ye_lo);
    ye_hi = amrex::max(zion[n] * aion_inv[n], ye_hi);
  }

  if (state.ye < ye_lo || state.ye > ye_hi){
    input_ye_is_valid = false;
  }
  
  if (!input_ye_is_valid) {
    // ensure Ye is valid
    composition(state);
  }

  // invoke newton-raphson to solve chemical potential of proton and neutron
  nse_nr_solver(state, eps);
  
  // get the nse_state
  T nse_state = get_nse_state(state);
  return nse_state;
    
}

#endif
