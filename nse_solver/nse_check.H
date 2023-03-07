#ifndef NSE_CHECK_H
#define NSE_CHECK_H

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
#include <nse_solver.H>

// Currently doesn't support aprox networks, only networks produced by pynucastro

#ifndef NEW_NETWORK_IMPLEMENTATION

// First check to see if we're in the ballpark of nse state

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void check_nse_molar(burn_t& state, const burn_t& nse_state, bool& nse_check){

  // This function gives the first estimate whether we're in the nse or not
  // it checks whether the molar fractions of n,p,a are approximately in NSE

  int count = 0;
  amrex::Real r = 1.0_rt;
  amrex::Real r_nse = 1.0_rt;

  // Check if neutron, proton, and helium are present in the network
  for (int n = 0; n < NumSpec; ++n){
    if (n == NSE_INDEX::p_index || n == NSE_INDEX::n_index|| n == NSE_INDEX::he4_index){
      ++count;
    }
  }

  if (count < 3){
    amrex::Error("Current network does not include the proton, neutron, or helium-4, thus cannot proceed with ASE.");
  }

  // Check if n,p,a are in equilibrium
  // these two ratios are defined in the ASE paper to determine whether network is in equilibrium

  for (int n = 0; n < NumSpec; ++n){

    if (n == NSE_INDEX::p_index || n == NSE_INDEX::n_index){
      r /= state.xn[n] * state.xn[n] * aion_inv[n] * aion_inv[n];
      r_nse /= nse_state.xn[n] * nse_state.xn[n] * aion_inv[n] * aion_inv[n];
    }

    else if (n == NSE_INDEX::he4_index){
      r *= state.xn[n] * aion_inv[n];
      r_nse *= nse_state.xn[n] * aion_inv[n];
    }
  }

  // equilibrium condition: if pass proceed with ase if not proceed with regular eos integration
  // Eq. 14 in Kushnir paper

  if (std::abs((r - r_nse) / r_nse) < 0.5_rt){
    nse_check = true;
  }
  else{
    nse_check = false;
  }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rate_by_nuc(int& rate_index, amrex::Array1D<int, 1, 3> reactants_indices,
		     amrex::Array1D<int, 1, 3> products_indices){

  rate_index = -1;
  bool found_index = false;

  // sort indices just to make sure

  sort_Array1D<1, 3>(reactants_indices);
  sort_Array1D<1, 3>(products_indices);

  // check for which rate matches these input indices
  // fill in the corresponding rate_index by providing array of reactants and products indices

  for (int i = 1; i <= Rates::NumRates; ++i){
    for (int j = 1; j <= 3; ++j){
      if (NSE_INDEX::rate_indices(i, j) == reactants_indices(j) && NSE_INDEX::rate_indices(i, j+3) == products_indices(j)){
	found_index = true;
      }
      else{
	found_index = false;
	break;
      }
    }

    if (found_index){
      rate_index = i;
      return;
    }
  }

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void find_fast_reaction_cycle(const int current_nuc_ind, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
			const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
			const burn_t& state, const amrex::Real& t_s, bool& found_fast_reaction_cycle){

  // This function finds the fast reaction cycle of a given nuclei.
  // This is done to satisfy Section 4.4.3 of Kushnir

  // Here pre-determine the max number of nuclei possible in each cycle
  // so we don't loop over too many nucs

  int max_nucs;
  find_max_nucs(max_nucs, current_nuc_ind);

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

  // initialize reactions and fast_reactions to -1

  for (int i = 1; i <= 6; ++i){
    for (int j = 1; j <= max_nucs; ++j){
      reactions(i, j) = -1;
      fast_reactions(i, j) = -1;
    }
  }

  // indices for reaction array to locate intermediate nuclei

  int inter_nuc_ind1 = 1;
  int inter_nuc_ind2 = 1;
  int inter_nuc_ind3 = 1;
  int inter_nuc_ind4 = 1;

  // the beginning nuclei of a reaction cycle is always itself

  reactions(1, 1) = current_nuc_ind;

  // Now store all possible nuclei involved in a reaction cycle started by nuclei: current_nuc_ind

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

  // Checks if there are intermediate nuclei in the cycle
  // if no intermediate nuclei present, then there is no way to make a full reaction cycle

  for (int i = 2; i <= 4; ++i){

    have_required_nucs = false;

    for (int j = 1; j <= max_nucs; ++j){
      if (reactions(i, j) != -1){
	have_required_nucs = true;
      }
    }

    if (!have_required_nucs){
      found_fast_reaction_cycle = false;

      return;
    }
  }


  // Now we have all the intermediate nuclei necessary for the reaction cycle
  // Now we check whether the possible reaction are fast enough

  // The first nuc in fast reaction is always itself

  fast_reactions(1, 1) = current_nuc_ind;
  amrex::Array1D<amrex::Real, 1, 2>  Y_i;

  // loop over reaction stages

  for (int i = 1; i <= 5; ++i){

    // loop over nuclei in the current stage

    for (int j = 1; j <= max_nucs; ++j){

      // skip this loop if current nuc is the same as the starting nuc
      // also skip if current index is -1, which is not valid

      if (reactions(i, j) == -1 || (reactions(i, j) == current_nuc_ind && i != 1)){
	continue;
      }

      // current current stage of nuc index, (i,j) to the next stage nuc index (i+1,k)
      // get nuc index of the next stage in reactions

      for (int k = 1; k <= max_nucs; ++k){

	// if reaction index is -1, skip, since its not a valid index

	if (reactions(i+1, k) == -1){
	  continue;
	}

	// determine which reaction it went through
	// Kushnir listed a total of 5 possible reactions to consider
	// (a,g), (a,p), (a,n), (g,p), (g,n)

	if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 4 &&
	    zion[reactions(i, j)] == zion[reactions(i+1, k)] - 2){

	  // (alpha, gamma) reaction
	  // implies reactants are he4, reactions(i, j) and product is (reactions(i+1, k)

	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = NSE_INDEX::he4_index;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = -1;
	  product_scratch(3) = -1;

	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);

	  // skip if there is no such rate.

	  if (rate_index == -1){
	    continue;
	  }

	  b_f = screened_rates(rate_index) * Y(reactions(i, j) + 1) * Y(NSE_INDEX::he4_index + 1) * state.rho;
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate

	  if (rate_index == -1){
	    continue;
	  }

	  b_r = screened_rates(rate_index) * Y(reactions(i+1, k) + 1);

	  Y_i(1) = Y(NSE_INDEX::he4_index + 1);
	  Y_i(2) = 0.0_rt;
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 3 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] - 2){

	  // (alpha, n) reaction

	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = NSE_INDEX::he4_index;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::n_index;
	  product_scratch(3) = -1;

	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);

	  if (rate_index == -1){
	    continue;
	  }

	  b_f = screened_rates(rate_index) * Y(reactions(i, j) + 1) * Y(NSE_INDEX::he4_index + 1) * state.rho;
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate

	  if (rate_index == -1){
	    continue;
	  }

	  b_r = screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::n_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::he4_index + 1);
	  Y_i(2) = Y(NSE_INDEX::n_index + 1);
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] - 3 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] - 1){

	  // (alpha, p) reaction

	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = NSE_INDEX::he4_index;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::p_index;
	  product_scratch(3) = -1;

	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);

	  if (rate_index == -1){
	    continue;
	  }

	  b_f = screened_rates(rate_index) * Y(reactions(i, j) + 1) * Y(NSE_INDEX::he4_index + 1) * state.rho;
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate

	  if (rate_index == -1){
	    continue;
	  }

	  b_r = screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::p_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::he4_index + 1);
	  Y_i(2) = Y(NSE_INDEX::p_index + 1);
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] + 1 &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)]){

	  // (gamma, n) reaction

	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = -1;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::n_index;
	  product_scratch(3) = -1;

	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);

	  if (rate_index == -1){
	    continue;
	  }

	  b_f = screened_rates(rate_index) * Y(reactions(i, j) + 1);
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate

	  if (rate_index == -1){
	    continue;
	  }

	  b_r = screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::n_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::n_index + 1);
	  Y_i(2) = 0.0_rt;
	}

	else if (aion[reactions(i, j)] == aion[reactions(i+1, k)] + 1  &&
		 zion[reactions(i, j)] == zion[reactions(i+1, k)] + 1){

	  // (gamma, p) reaction

	  reaction_scratch(1) = reactions(i, j);
	  reaction_scratch(2) = -1;
	  reaction_scratch(3) = -1;

	  product_scratch(1) = reactions(i+1, k);
	  product_scratch(2) = NSE_INDEX::p_index;
	  product_scratch(3) = -1;

	  get_rate_by_nuc(rate_index, reaction_scratch, product_scratch);

	  if (rate_index == -1){
	    continue;
	  }

	  b_f = screened_rates(rate_index) * Y(reactions(i, j) + 1);
	  get_rate_by_nuc(rate_index, product_scratch, reaction_scratch);

	  // skip if there are no both forward and reverse rate

	  if (rate_index == -1){
	    continue;
	  }

	  b_r = screened_rates(rate_index) * Y(reactions(i+1, k) + 1) * Y(NSE_INDEX::p_index + 1) * state.rho;

	  Y_i(1) = Y(NSE_INDEX::p_index + 1);
	  Y_i(2) = 0.0_rt;
	}

	else{

	  // continue if no matching reactions

	  continue;
	}

	// A condition to check whether rates are fast enough given by Eq. 11 in Kushnir.
	// Not sure if the 3rd condition (Eq. 12) is required
	// It didn't mention it in in the paragraph in section 4.4.3, but listed it in Eq. 12
	// It didn't have a significant difference, but I'll leave it here for future reference

	if (
	    ((2.0_rt * std::abs(b_f - b_r) / (b_f + b_r)) < ase_tol)
	    &&
	    ((Y_i(1) / amrex::min(b_f, b_r) < ase_tol * t_s)
	     || (Y_i(2) / amrex::min(b_f, b_r) < ase_tol * t_s))
	  ) {

	  // store nuclei index to fast_reaction if pass the condition

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

	// if we see non -1 index, which means we found nucs in the next stage
	// that are fast enough, but not yet complete the full cycle

	if (fast_reactions(i + 1, p) != -1){
	  have_required_nucs = true;
	}
    }

    // if there is not a single nuc in the next stage meet the requirement
    // then there is no intermdiate nuclei, then break.

    if (!have_required_nucs){
      found_fast_reaction_cycle = false;
      return;
    }
  }

}


// After all preliminary checks are successful, lets do nse grouping.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
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


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_union(const int nuc_ind_a, const int nuc_ind_b, amrex::Array1D<int, 1, NumSpec>& group_ind){

  // This function joins the two group of the two nuc indices:nuc_ind_a and nuc_ind_b
  // The smaller group is joined to the larger group.

  int root_index_a = get_root_index(nuc_ind_a, group_ind);
  int root_index_b = get_root_index(nuc_ind_b, group_ind);

  if (root_index_a == root_index_b){
    return;
  }

  // find size of the two groups containing a and b

  int group_a_size = 0;
  int group_b_size = 0;

  for (int n = 0; n < NumSpec; ++n){

    if (get_root_index(n, group_ind) == root_index_a){
      ++group_a_size;
    }
    else if (get_root_index(n, group_ind) == root_index_b){
      ++group_b_size;
    }
  }

  // merge group with less isotopes to group with larger isotope

  if (group_a_size >= group_b_size){
    group_ind(root_index_b) = group_ind(root_index_a);
  }

  else{
    group_ind(root_index_a) = group_ind(root_index_b);
  }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
int get_pair_rate_index(const int& rate_index){

  // This funciton finds the pair rate index given a rate.
  // For example given forward rate index, return reverse rate index, vice versa

  int pair_rate_index = -1;
  bool found_reverse_rate = false;
  bool valid_rate = false;

  // if all -1 then return, since its an invalid rate

  for (int n = 1; n <= 6; ++n){
    if (NSE_INDEX::rate_indices(rate_index, n) != -1){
      valid_rate = true;
      break;
    }
  }

  if (!valid_rate){
    return pair_rate_index;
  }


  for (int n = 1; n <= Rates::NumRates; ++n) {
    for (int i = 1; i <= 3; ++i) {
      if (NSE_INDEX::rate_indices(n, i) == NSE_INDEX::rate_indices(rate_index, i+3)
	  && NSE_INDEX::rate_indices(n, i+3) == NSE_INDEX::rate_indices(rate_index, i)
	  ) {
	found_reverse_rate = true;
      }
      else {
	found_reverse_rate = false;
	break;
      }
    }

    if (found_reverse_rate) {
      pair_rate_index = n;
      break;
    }
  }

  return pair_rate_index;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_single_group(const amrex::Array1D<int, 1, NumSpec>& group_ind){

  // This function checks whether all isotopes are either in the LIG group
  // or in another single group.

  int nonLIG_index = -1;
  int LIG_root_index = get_root_index(NSE_INDEX::he4_index, group_ind);

  for (int n = 0; n < NumSpec; ++n){

    if (get_root_index(n, group_ind) == LIG_root_index){
      continue;
    }

    if (nonLIG_index == -1){
      nonLIG_index = get_root_index(n, group_ind);
    }

    else if (get_root_index(n, group_ind) != nonLIG_index){
      return false;
    }

  }

  return true;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void is_fastest_rate(amrex::Array1D<int, 1, 2>& merge_indices,
		     amrex::Real& fastest_t, const int current_rate,
		     const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, const burn_t& state,
		     const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
		     const amrex::Array1D<int, 1, NumSpec>& group_ind,
		     const amrex::Real& t_s) {

  // This function determines whether current_rate has a faster timescale than
  // the provided time scale, current fastest_t.
  // Updates the timescale once if the current rate has a faster timescale

  // It also skips current_rate if certain condition is met.

  int root_index;
  int num_nonLIG = 0;
  int nonLIG_root = -1;

  // Find the Y_group of the nuclei involved in the net
  // Y_group are the group molar fractions corresponding to the reactants and products

  // first 3 is for reactants, last 3 for products

  amrex::Array1D<amrex::Real, 1, 6> Y_group = {0.0_rt, 0.0_rt, 0.0_rt,
					       0.0_rt, 0.0_rt, 0.0_rt};

  // Here Y_nonNPA represent 2 group isotopes that are not n,p,a
  // there can be at most 2 nonNPA isotopes in a reaction for it to be valid

  amrex::Array1D<int, 1, 2> nonNPA_ind = {-1, -1};
  int m = 1;                    // m is counts the number of nonNPA isotope detected.

  // find the reverse rate index of the current rate

  int pair_rate_index = get_pair_rate_index(current_rate);

  // Check if there are > 2 isotopes or 0 isotopes in LIG (light isotope group) in this reaction

  for (int k = 1; k <= 6; ++k){

    // skip if current isotope is -1, which means null

    if (NSE_INDEX::rate_indices(current_rate, k) == -1){
      continue;
    }

    // Find the root index of this isotope

    root_index = get_root_index(NSE_INDEX::rate_indices(current_rate, k), group_ind);

    // Determine number of nonLIG isotopes
    // also check whether nonLIG isotopes are already merged

    if (root_index != get_root_index(NSE_INDEX::he4_index, group_ind)){

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

    // let LIG group use their own Y instead of Y_group. (not sure if this is true)

    if (root_index == get_root_index(NSE_INDEX::he4_index, group_ind)){

	Y_group(k) = Y(NSE_INDEX::rate_indices(current_rate, k) + 1);
    }

    else{
      for (int n = 0; n < NumSpec; ++n){

	// if nuclei have the same root_index, then they are the same group

	if (get_root_index(n, group_ind) == root_index){
	  Y_group(k) += Y(n + 1);
	}

      }
    }

    // Find indices that are non n,p,a

    if (NSE_INDEX::rate_indices(current_rate, k) != NSE_INDEX::p_index
	&& NSE_INDEX::rate_indices(current_rate, k) != NSE_INDEX::n_index
	&& NSE_INDEX::rate_indices(current_rate, k) != NSE_INDEX::he4_index
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
  // 2) skip when reaction has more than 3 reactants or products involved
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

  b_f = screened_rates(current_rate) * Y(NSE_INDEX::rate_indices(current_rate, 3) + 1);
  b_r = screened_rates(pair_rate_index) * Y(NSE_INDEX::rate_indices(current_rate, 6) + 1);

  if (NSE_INDEX::rate_indices(current_rate, 2) != -1){

    if (NSE_INDEX::rate_indices(current_rate, 2) == NSE_INDEX::rate_indices(current_rate, 3)){
      b_f *= 0.5_rt;
    }

    b_f *= Y(NSE_INDEX::rate_indices(current_rate, 2) + 1) * state.rho;
  }

  if (NSE_INDEX::rate_indices(current_rate, 5) != -1){

    if (NSE_INDEX::rate_indices(current_rate, 5) == NSE_INDEX::rate_indices(current_rate, 6)){
      b_r *= 0.5_rt;
    }

    b_r *= Y(NSE_INDEX::rate_indices(current_rate, 5) + 1) * state.rho;
  }

  // Find the timescale of the rate, See Eq. 17 and Eq. 11 in Kushnir

  amrex::Real scratch_t = Y_group(nonNPA_ind(1)) / amrex::min(b_f, b_r);
  if (nonNPA_ind(2) != -1){
    scratch_t = amrex::min(scratch_t, Y_group(nonNPA_ind(2)) / amrex::min(b_f, b_r));
  }

  // return if current rate timescale is larger than previous time scale

  if (fastest_t < scratch_t){
    return;
  }

  // if current rate timescale is smaller (faster), have additional checks
  // to see whether the forward and balance rates are in equilibrium
  // see Eq. 17 again

  if ( (2.0_rt * std::abs(b_f - b_r) / (b_f + b_r) < ase_tol)
      && (scratch_t < ase_tol * t_s) ) {

    fastest_t = scratch_t;

    // get merge indices for nonNPA isotopes

    for (int n = 1; n <= 2; ++n){

      // If nonNPA index is -1, meaning null, then use LIG index

      if (nonNPA_ind(n) == -1){
	merge_indices(n) = get_root_index(NSE_INDEX::he4_index, group_ind);
      }
      else{
	merge_indices(n) = NSE_INDEX::rate_indices(current_rate, nonNPA_ind(n));
      }
    }

  }

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_grouping(amrex::Array1D<int, 1, NumSpec>& group_ind, const burn_t& state,
		  const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
		  const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
		  const amrex::Real& t_s){

  // This function groups all the nuclei using group_ind
  // which contains the node #

  // fill in initial group_ind, group_ind go from 1 to NumSpec

  for (int i = 1; i <= NumSpec; ++i){
    group_ind(i) = i;

    // let n,p,a form the same group (LIG) initially, let 1 be index of LIG

    if (i == NSE_INDEX::p_index + 1 || i == NSE_INDEX::n_index + 1
	|| i == NSE_INDEX::he4_index + 1){
      group_ind(i) = 1;
    }
  }

  amrex::Array1D<int, 1, 2> merge_indices;
  amrex::Real fastest_t;

  // Go through each reaction and perform grouping

  while(true){

    merge_indices(1) = -1;
    merge_indices(2) = -1;
    fastest_t = std::numeric_limits<amrex::Real>::max();

    // Find the fastest time scale of the avaliable reaction

    for (int n = 1; n <= Rates::NumRates; ++n){

      // Check if current rate is the faster rate than previous rate

      is_fastest_rate(merge_indices, fastest_t, n, Y, state, screened_rates, group_ind, t_s);
    }

    // if there are no reaction satisfy conditions, break

    if (merge_indices(1) == -1 && merge_indices(2) == -1){
      break;
    }

    // union the isotopes of the rate that have the fastest time scale

    nse_union(merge_indices(1), merge_indices(2), group_ind);
  }

}
#endif


AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_nse(burn_t& current_state, bool skip_molar_check=false) {

  // This function returns the boolean that tells whether we're in nse or not
  // Note that it only works with pynucastro network for now.
  
#ifndef NEW_NETWORK_IMPLEMENTATION

  current_state.nse = false;

  if (current_state.T < 2.0e9_rt) {
    return current_state.nse;
  }
    
  // Get the nse state which is used to compare nse molar fractions.

  const auto nse_state = get_actual_nse_state(current_state);

  // Check whether state is in the ballpark of NSE

  if (!skip_molar_check) {

    check_nse_molar(current_state, nse_state, current_state.nse);    
    if (!current_state.nse) {
      return current_state.nse;
    }
  }

  // We can do a further approximation where we use the NSE mass fractions
  // instead of the current mass fractions. This makes the check soley dependent on
  // the thermodynamic condition.
  // Note we only do this after the first check, which should tell us whether
  // our current mass fractions are in the ballpark of NSE mass fractions.

  burn_t state = current_state;

  if (nse_molar_independent) {
    state = nse_state;
    state.dx = current_state.dx;
  }

  // set molar fractions
  // initialize ydot to store Ydot and energy generation rate.

  amrex::Array1D<amrex::Real, 1, NumSpec> Y;
  amrex::Array1D<amrex::Real, 1, NumSpec + 1> ydot;

  for (int n = 1; n <= NumSpec; ++n){
    Y(n) = state.xn[n-1] * aion_inv[n-1];
    ydot(n) = 0.0_rt;
  }

  rate_t rate_eval;
  constexpr int do_T_derivatives = 0;
  evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

  eos_t eos_state;                             // need eos_state for speed of sound

  // Initialize t_s, sound crossing timescale for a single zone to be max

  amrex::Real t_s = std::numeric_limits<amrex::Real>::max();

  // If we care about checking the timescale of the rate to be smaller than t_s, then:

  if (!nse_dx_independent) {
    burn_to_eos<eos_t>(state, eos_state);
    eos(eos_input_rt, eos_state);
    t_s = state.dx / eos_state.cs;   // a parameter to characterize whether a rate is fast enough
  }

  // Find ydot

  rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

  // Find energy generation rate

  amrex::Real enuc;
  ener_gener_rate(ydot, enuc);

#ifdef NEUTRINOS
  // get abar and zbar
  composition(state);
  
  amrex::Real sneut, dsneutdt, dsneutdd, snuda, snudz;
  sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);
#else
  amrex::Real sneut = 0.0_rt, dsneutdt =0.0_rt , dsneutdd = 0.0_rt, snuda = 0.0_rt, snudz = 0.0_rt;
#endif

  // fill in energy generation rate to ydot

  ydot(NumSpec + 1) = enuc - sneut;

  // Now we look through the network and see if there are fast reaction cycles
  // Need to separate forward and reverse rate and detemine each step is fast enough.
  // use vectors for now

  bool found_fast_reaction_cycle = false;

  // Do a reverse for loop to start from heaviest nuclei

  for (int n = NumSpec-1; n >= 0; --n){
    if (found_fast_reaction_cycle){
      break;
    }

    find_fast_reaction_cycle(n, Y, rate_eval.screened_rates, state, t_s, found_fast_reaction_cycle);
  }

  // Do nse grouping if found fast reaction cycle.

  // This holds group index for each nuclei
  // it determines which group each nuclei is in

  amrex::Array1D<int, 1, NumSpec> group_ind;

  // if not found_fast_reaction_cycle, state not in nse

  if (!found_fast_reaction_cycle){
    current_state.nse = false;
  }
  else{

    // Do nse grouping if found_fast_reaction_cycle

    nse_grouping(group_ind, state, Y, rate_eval.screened_rates, t_s);

    // Check if we result in a single group after grouping

    current_state.nse = false;

    if (in_single_group(group_ind)){
      current_state.nse = true;
    }
  }

  return current_state.nse;

#else
  amrex::Error("in_nse() is currently not supported for aprox networks!");

  return false;
#endif

}

#endif
