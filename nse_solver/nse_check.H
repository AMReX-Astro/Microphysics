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
#include <nse_solver.H>

// Currently doesn't support aprox networks, only networks produced by pynucastro

#ifndef NEW_NETWORK_IMPLEMENTATION

// First check to see if we're in the ballpark of nse state

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void check_nse_molar(const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                     const amrex::Array1D<amrex::Real, 1, NumSpec>& Y_nse,
                     bool& nse_check) {

    // This function gives the first estimate whether we're in the nse or not
    // it checks whether the molar fractions of n,p,a are approximately in NSE

    amrex::Real r = 1.0_rt;
    amrex::Real r_nse = 1.0_rt;

    nse_check = false;

    // raise error if no proton or helium-4 in the network.

    if (NSE_INDEX::H1_index == -1 || NSE_INDEX::He4_index == -1) {
        amrex::Error("Need proton and helium-4 in the network for NSE_NET to work");
    }

    // If there are proton, or helium in the network

    // Check if n,p,a are in equilibrium
    // these two ratios are defined in the ASE paper to determine whether network is in equilibrium

    for (int n = 0; n < NumSpec; ++n) {

        if (n == NSE_INDEX::H1_index || n == NSE_INDEX::N_index) {
            r /= Y(n+1) * Y(n+1);
            r_nse /= Y_nse(n+1) * Y_nse(n+1);
        }

        else if (n == NSE_INDEX::He4_index) {
            r *= Y(n+1);
            r_nse *= Y_nse(n+1);
        }
    }

    // equilibrium condition: if pass proceed with ase if not proceed with regular eos integration
    // Eq. 14 in Kushnir paper

    // if there is neutron in the network

    if ((std::abs(r - r_nse) < 0.5_rt * r_nse)
        && (NSE_INDEX::N_index != -1)) {
        nse_check = true;
        return;
    }

    // if there is no neutron in the network

    if ((std::abs(r - r_nse) < 0.25_rt * r_nse)
        && (NSE_INDEX::N_index == -1)) {
        nse_check = true;
        return;
    }

    // Overall molar fraction check

    for (int n = 0; n < NumSpec; ++n) {
        amrex::Real abs_diff = std::abs(Y(n+1) - Y_nse(n+1));
        if (abs_diff > nse_abs_tol && abs_diff > nse_rel_tol * Y(n+1)) {
            return;
        }
    }
    nse_check = true;
}


// After all preliminary checks are successful, lets do nse grouping.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
int get_root_index(const int nuc_ind,
                   const amrex::Array1D<int, 1, NumSpec>& group_ind) {

    // This function returns the root index of the nuclei
    // by providing the nuclei index [0, NumSpec-1], and group indices, group_ind

    int root_index;
    int scratch_ind = nuc_ind;

    while(true) {
        root_index = group_ind(scratch_ind + 1);

        if (root_index != scratch_ind + 1) {
            scratch_ind = root_index - 1;
        }

        else {
            return root_index;
        }
    }

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_union(const int nuc_ind_a, const int nuc_ind_b, amrex::Array1D<int, 1, NumSpec>& group_ind) {

    // This function joins the two group of the two nuc indices:nuc_ind_a and nuc_ind_b
    // The smaller group is joined to the larger group.

    int root_index_a = get_root_index(nuc_ind_a, group_ind);
    int root_index_b = get_root_index(nuc_ind_b, group_ind);

    if (root_index_a == root_index_b) {
        return;
    }

    // find size of the two groups containing a and b

    int group_a_size = 0;
    int group_b_size = 0;

    for (int n = 0; n < NumSpec; ++n) {

        if (get_root_index(n, group_ind) == root_index_a) {
            ++group_a_size;
        }
        else if (get_root_index(n, group_ind) == root_index_b) {
            ++group_b_size;
        }
    }

    // merge group with less isotopes to group with larger isotope

    if (group_a_size >= group_b_size) {
        group_ind(root_index_b) = group_ind(root_index_a);
    }
    else {
        group_ind(root_index_a) = group_ind(root_index_b);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_single_group(const amrex::Array1D<int, 1, NumSpec>& group_ind) {

    // This function checks whether all isotopes are either in the LIG group
    // or in another single group.

    int LIG_root_index = get_root_index(NSE_INDEX::He4_index, group_ind);

    int nonLIG_index = -1;

    int oddN_group = -1;
    int evenN_group = -1;

    bool in_single_group = true;

    // Consider NSE when there is a single group with an optional LIG group

    for (int n = 0; n < NumSpec; ++n) {
        if (get_root_index(n, group_ind) == LIG_root_index) {
            continue;
        }

        if (nonLIG_index == -1) {
            nonLIG_index = get_root_index(n, group_ind);
            continue;
        }

        if (get_root_index(n, group_ind) != nonLIG_index) {
            in_single_group = false;
            break;
        }

    }

    // If there no neutrons are in the network and original condition failed
    // Consider a looser condition by looking at nuclei heavier than Si28
    // There seems to be two big groups after Si28 in NSE:
    // 1) isotopes with even N
    // 2) isotopes with odd N

    if (NSE_INDEX::N_index == -1 && !in_single_group) {

        in_single_group = true;

        for (int n = 0; n < NumSpec; ++n) {

            if (zion[n] >= 14) {

                // Get even N group index
                if (evenN_group == -1 && std::fmod(aion[n] - zion[n], 2) == 0.0_rt) {
                    evenN_group = get_root_index(n, group_ind);
                    continue;
                }

                // Get odd N group index
                if (oddN_group == -1 && std::fmod(aion[n] - zion[n], 2) == 1.0_rt) {
                    oddN_group = get_root_index(n, group_ind);
                    continue;
                }

                if ((std::fmod(aion[n] - zion[n], 2) == 0.0_rt && evenN_group != get_root_index(n, group_ind)) ||
                    (std::fmod(aion[n] - zion[n], 2) == 1.0_rt && oddN_group != get_root_index(n, group_ind))) {
                    in_single_group = false;
                    break;
                }
            }

        }
    }

    return in_single_group;
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fill_reaction_timescale(amrex::Array1D<T, 1, Rates::NumRates>& reaction_timescales,
                             const int current_rate_index, const amrex::Real rho,
                             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                             const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
                             const amrex::Real t_s) {
    // This function fills in the timescale of the reaction for the rates that
    // are compatible for NSE.

    // Default to the largest possible timescale

    constexpr amrex::Real max_timescale = std::numeric_limits<amrex::Real>::max();
    reaction_timescales(current_rate_index) = max_timescale;

    //
    // Few conditions to check:
    // 1) skip if there is no reverse rate involved.
    // 2) skip when reaction has more than 3 reactants or products involved
    // If these conditions are met, we don't consider them for NSE
    // so we assume they have the largest (slowest) timescale
    //

    int reverse_rate_index = NSE_INDEX::rate_indices(current_rate_index, 7);
    if ( (reverse_rate_index == -1)
         || (NSE_INDEX::rate_indices(current_rate_index, 1) != -1)
         || (NSE_INDEX::rate_indices(current_rate_index, 4) != -1)
        ) {
        return;
    }

    //
    // 3) If there are more than 2 non Neutron, Proton, or Helium-4 in the rate
    //

    amrex::Array1D<int, 1, 2> non_NHA_ind = {-1, -1};
    int non_NHA_counts = 0;
    for (int k = 2; k <= 6; ++k) {
        if (NSE_INDEX::rate_indices(current_rate_index, k) == -1) {
            continue;
        }

        int is_neutron_in_network = (NSE_INDEX::N_index != -1);
        int is_valid_index = (NSE_INDEX::rate_indices(current_rate_index, k) != NSE_INDEX::H1_index &&
                              NSE_INDEX::rate_indices(current_rate_index, k) != NSE_INDEX::He4_index);

        if ((!is_neutron_in_network && is_valid_index) ||
            (is_neutron_in_network && is_valid_index &&
             NSE_INDEX::rate_indices(current_rate_index, k) != NSE_INDEX::N_index)) {

            ++non_NHA_counts;

            // Check if count exceeds 2

            if (non_NHA_counts > 2) {
                return;
            }

            // Store the index where we have non Neutron, Proton, or Helium-4

            non_NHA_ind(non_NHA_counts) = k;
        }
    }

    // Calculate the forward and reverse rates of the current rate index

    amrex::Real b_f;
    amrex::Real b_r;

    b_f = screened_rates(current_rate_index) * Y(NSE_INDEX::rate_indices(current_rate_index, 3) + 1);
    b_r = screened_rates(reverse_rate_index) * Y(NSE_INDEX::rate_indices(current_rate_index, 6) + 1);

    if (NSE_INDEX::rate_indices(current_rate_index, 2) != -1) {
        if (NSE_INDEX::rate_indices(current_rate_index, 2) == NSE_INDEX::rate_indices(current_rate_index, 3)) {
            b_f *= 0.5_rt;
        }
        b_f *= Y(NSE_INDEX::rate_indices(current_rate_index, 2) + 1) * rho;
    }

    if (NSE_INDEX::rate_indices(current_rate_index, 5) != -1) {
        if (NSE_INDEX::rate_indices(current_rate_index, 5) == NSE_INDEX::rate_indices(current_rate_index, 6)) {
            b_r *= 0.5_rt;
        }
        b_r *= Y(NSE_INDEX::rate_indices(current_rate_index, 5) + 1) * rho;
    }

    // Find the timescale of the rate, See Equation 11 in Kushnir

    //
    // Note that here I made a simplification compared to Kushnir:
    // When they calculate the reaction timescale,
    // they used the molar fraction of the entire group, but I just
    // use the molar fraction of isotopes participated in the rate itself.
    // This complicates calculation since Y_group changes as you merge
    // isotopes into groups, so we need to do calculations iteratively.
    //

    amrex::Real t_i = Y(NSE_INDEX::rate_indices(current_rate_index, non_NHA_ind(1))+1) / amrex::min(b_f, b_r);
    if (non_NHA_ind(2) != -1) {
        t_i = amrex::min(t_i, Y(NSE_INDEX::rate_indices(current_rate_index, non_NHA_ind(2))+1) / amrex::min(b_f, b_r));
    }

    //
    // Condition for checking if forward and reverse rates are in equilibrium
    // Also whether the reaction time scale is smaller than sound crossing time
    // See Equation 17 in Kushnir Paper
    //

    if ((2.0_rt * std::abs(b_f - b_r) / (b_f + b_r) < ase_tol) &&
        (t_i < ase_tol * t_s)) {
        reaction_timescales(current_rate_index) = t_i;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fill_merge_indices(amrex::Array1D<int, 1, 2>& merge_indices,
                        const int current_rate_index,
                        const amrex::Array1D<int, 1, NumSpec>& group_ind) {
    // This function determines the merge indices for the current rate index

    // First determine the non neutron, proton, helium-4 indices of the rate

    amrex::Array1D<int, 1, 2> non_NHA_ind = {-1, -1};
    int non_NHA_counts = 0;

    merge_indices(1) = -1;
    merge_indices(2) = -1;

    for (int k = 2; k <= 6; ++k) {
        if (NSE_INDEX::rate_indices(current_rate_index, k) == -1) {
            continue;
        }

        int is_neutron_in_network = (NSE_INDEX::N_index != -1);
        int is_valid_index = (NSE_INDEX::rate_indices(current_rate_index, k) != NSE_INDEX::H1_index &&
                              NSE_INDEX::rate_indices(current_rate_index, k) != NSE_INDEX::He4_index);

        if ((!is_neutron_in_network && is_valid_index) ||
            (is_neutron_in_network && is_valid_index &&
             NSE_INDEX::rate_indices(current_rate_index, k) != NSE_INDEX::N_index)) {

            ++non_NHA_counts;

            if (non_NHA_counts > 2) {
                return;
            }

            // Store the index where we have non Neutron, Proton, or Helium-4

            non_NHA_ind(non_NHA_counts) = k;
        }
    }

    // Check whether isotopes are already merged

    int num_nonLIG = 0;
    int nonLIG_root = -1;

    for (int k = 2; k <= 6; ++k) {
        if (NSE_INDEX::rate_indices(current_rate_index, k) == -1) {
            continue;
        }

        int root_index = get_root_index(NSE_INDEX::rate_indices(current_rate_index, k),
                                        group_ind);

        // Determine number of nonLIG isotopes
        // also check whether nonLIG isotopes are already merged

        if (root_index != get_root_index(NSE_INDEX::He4_index, group_ind)) {

            ++num_nonLIG;

            //
            // return if nonLIG_root index is repeated, i.e. isotopes already merged
            // Initialize nonLIG_root = -1.
            // Update nonLIG_root during first encounter nonLIG
            // If nonLIG_root is repeated during second encounter, then abort
            //

            if (root_index == nonLIG_root) {
                return;
            }

            // Update the nonLIG_root during first encounter

            nonLIG_root = root_index;
        }
    }

    // skip if number of LIG is greater than 2 or equal to 0

    if (num_nonLIG == 0 || num_nonLIG > 2) {
        return;
    }

    for (int n = 1; n <= 2; ++n) {

        // If non_NHA index is -1, meaning null, then use LIG index

        if (non_NHA_ind(n) == -1) {
            merge_indices(n) = get_root_index(NSE_INDEX::He4_index, group_ind);
        }
        else {
            merge_indices(n) = NSE_INDEX::rate_indices(current_rate_index, non_NHA_ind(n));
        }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_grouping(amrex::Array1D<int, 1, NumSpec>& group_ind, const amrex::Real rho,
                  const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                  const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
                  const amrex::Real t_s) {

    // This function groups all the nuclei using group_ind
    // which contains the node #

    // fill in initial group_ind, group_ind go from 1 to NumSpec

    for (int i = 1; i <= NumSpec; ++i) {
        group_ind(i) = i;

        // let n,p,a form the same group (LIG) initially, let 1 be index of LIG

        if (i == NSE_INDEX::H1_index + 1 || i == NSE_INDEX::N_index + 1
            || i == NSE_INDEX::He4_index + 1) {

            // group_ind(i) = NSE_INDEX::He4_index;

            group_ind(i) = 1;
        }
    }

    // Let's first create an array for reaction_timescales and reaction_indices
    // Then fill in the reaction timescale and index for each rate.

    amrex::Array1D<amrex::Real, 1, Rates::NumRates> reaction_timescales;
    amrex::Array1D<int, 1, Rates::NumRates> rate_indices;

    for (int n = 1; n <= Rates::NumRates; ++n) {
        fill_reaction_timescale(reaction_timescales, n, rho, Y,
                                screened_rates, t_s);
        rate_indices(n) = n;
    }

    //
    // Sort rate_indices using reaction_timescales.
    // The sorted rate_indices should correspond to reactions
    // from smallest (fastest) to largest (slowest) timescale
    //

    quickSort_Array1D(rate_indices, reaction_timescales);

    // After the rate indices are sorted based on reaction timescales.
    // Now do the grouping based on the timescale.

    amrex::Array1D<int, 1, 2> merge_indices;
    constexpr amrex::Real max_timescale = std::numeric_limits<amrex::Real>::max();

    for (int n = 1; n <= Rates::NumRates; ++n) {

        int current_rate_index = rate_indices(n);

        //
        // Check if the timescale is at max, which means this reaction rate is
        // not valid for considering merging.
        //

        if (reaction_timescales(n) == max_timescale) {
            break;
        }

        // Fill in the merge index

        fill_merge_indices(merge_indices, current_rate_index, group_ind);

        // Check if merge indices are -1, which means nothing to merge

        if (merge_indices(1) == -1 || merge_indices(2) == -1) {
            continue;
        }

        // union the isotopes into the same group

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

    amrex::Real T_in = current_state.T_fixed > 0.0_rt ? current_state.T_fixed : current_state.T;

    // If temperature is below T_min_nse
    // Or if we enable direct by a simple temperature threshold,
    // so T_nse_net > 0.0 but T_in < T_nse_net, we abort early to avoid computing cost.

    if ((T_in < T_min_nse) || (T_nse_net > 0.0_rt && T_in < T_nse_net)) {
        return current_state.nse;
    }

    // Get the NSE state using (rho, T, Ye) as input
    // Then compare input molar fractions to the NSE molar fractions.

    const auto nse_state = get_actual_nse_state(nse_input_rty, current_state);

    // Convert to molar fractions

    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    amrex::Array1D<amrex::Real, 1, NumSpec> Y_nse;

    for (int n = 0; n < NumSpec; ++n) {
        Y(n+1) = current_state.y[SFS+n] * aion_inv[n] / current_state.rho;
        Y_nse(n+1) = nse_state.xn[n] * aion_inv[n];
    }

    // Check whether state is in the ballpark of NSE

    if (!skip_molar_check) {
        check_nse_molar(Y, Y_nse, current_state.nse);
        if (!current_state.nse) {
            return current_state.nse;
        }
    }

    // A simple temperature criteria after molar fraction check for determining NSE state
    // By default, T_nse_net = -1.0
    // So this is only enabled if the user provides value in the input file

    if (T_nse_net > 0.0_rt && T_in > T_nse_net) {
        current_state.nse = true;
        return current_state.nse;
    }

    // We can do a further approximation where we use the NSE mass fractions
    // instead of the current mass fractions. This makes the check solely dependent on
    // the thermodynamic condition.
    // Note we only do this after the first check, which should tell us whether
    // our current mass fractions are in the ballpark of NSE mass fractions.

    if (nse_molar_independent || skip_molar_check) {
        for (int n = 1; n <= NumSpec; ++n) {
            Y(n) = Y_nse(n);
        }
    }

    auto state = current_state;

    // Find the mass fraction of the current state

    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = Y(n+1) * aion[n];
    }

    rate_t rate_eval;
    constexpr int do_T_derivatives = 0;
    evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

    // need eos_state for speed of sound

    eos_t eos_state;

    // Initialize t_s, sound crossing timescale for a single zone to be max

    amrex::Real t_s = std::numeric_limits<amrex::Real>::max();

    // If we care about checking the timescale of the rate to be smaller than t_s, then:

    if (!nse_dx_independent) {
        burn_to_eos(state, eos_state);
        eos(eos_input_rt, eos_state);

        // a parameter to characterize whether a rate is fast enough

        t_s = state.dx / eos_state.cs;
    }

    amrex::Array1D<int, 1, NumSpec> group_ind;

    //
    // In Kushnir paper, we also need to perform a check
    // on finding at least 1 "fast reaction cycle"
    // However, it only works with neutron in the network.
    // I've also skipped this step when working with non-neutron
    // network. And it all seem to work fine.
    //

    // Now do nse grouping

    nse_grouping(group_ind, state.rho, Y, rate_eval.screened_rates, t_s);

    // Check if we result in a single group after grouping

    current_state.nse = in_single_group(group_ind);
    return current_state.nse;

#else
    amrex::Error("in_nse() is currently not supported for aprox networks!");
    return false;
#endif

}

#endif
