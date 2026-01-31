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


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void check_nse_molar(const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                     const amrex::Array1D<amrex::Real, 1, NumSpec>& Y_nse,
                     bool& nse_check) {

    // This function gives the first estimate whether we're in the nse or not
    // it checks whether the molar fractions of n,p,a are approximately in NSE
    // The idea is to consider a net reaction of the form: 2n + 2 H1 ⇄ He4
    //
    // The equilibrium condition enforces this reaction rate to take condition:
    // Y_{He4} / (Y_{n}^2 * Y_{H1}^2) = Y_{He4,NSE} / (Y_{n,NSE}^2 * Y_{H1,NSE}^2)
    // See Eq. 14 from Kushnir & Katz, 2019.
    //
    // If ratio check failed, then do an overall molar fraction check

    nse_check = false;

    if constexpr (NSE_INDEX::He4_index == -1) {
        amrex::Error("NSE Check requires at least He4 to be in the network.");
    }

    if constexpr (NSE_INDEX::N_index != -1 && NSE_INDEX::H1_index != -1) {

        // Molar fraction check from Eq. 14 in Kushnir & Katz, 2019.
        constexpr auto n = NSE_INDEX::N_index;
        constexpr auto p = NSE_INDEX::H1_index;
        constexpr auto a = NSE_INDEX::He4_index;

        if (Y(n) != 0.0_rt && Y(p) != 0.0_rt) {
            amrex::Real r = Y(a) / (Y(p) * Y(p) * Y(n) * Y(n));
            amrex::Real r_nse = Y_nse(a) / (Y_nse(p) * Y_nse(p) * Y_nse(n) * Y_nse(n));

            if (std::abs(r - r_nse) < 0.5_rt * r_nse) {
                nse_check = true;
                return;
            }
        }

    }

    // Do overall molar fraction check if the ratio check failed.

    for (int n = 1; n <= NumSpec; ++n) {
        if (std::abs(Y(n) - Y_nse(n)) > nse_rel_tol * Y_nse(n)) {
            return;
        }
    }
    nse_check = true;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
int get_root_index(const int nuc_idx, amrex::Array1D<int, 1, NumSpec>& group_idx) {
    // Returns the root index of the nuclei given the nuclei index [1, NumSpec]
    // Also updates all parent indices to point to the root index in group_idx
    //

    int scratch = nuc_idx;
    int root = scratch;

    // Find the root index.
    // Keep looping until the index of the current element
    // is the same as the value stored -- i.e. root element.
    while (group_idx(root) != root) {
        root = group_idx(root);
    }

    // Path compression.
    // Updates the group_idx so that all the parent nodes
    // will store the value that points directly to root.
    // This can make future get_root_index() faster.
    while (scratch != root) {
        int parent = group_idx(scratch);
        group_idx(scratch) = root;
        scratch = parent;
    }

    return root;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_union(const int nuc_idx_a, const int nuc_idx_b,
               amrex::Array1D<int, 1, NumSpec>& group_idx,
               amrex::Array1D<int, 1, NumSpec>& group_size) {
    // Joins two group represented by two nuclei index: nuc_idx_a and nuc_idx_b
    // Note that the smaller group is merged to the larger group.
    // group_size of the larger group is updated.

    int root_idx_a = get_root_index(nuc_idx_a, group_idx);
    int root_idx_b = get_root_index(nuc_idx_b, group_idx);

    // Return if they're the same group
    if (root_idx_a == root_idx_b) {
        return;
    }

    // Merge smaller group to the larger group.
    // Do the union and update group size of the larger group.
    // Set group_size of the merged group to be 0
    if (group_size(root_idx_a) >= group_size(root_idx_b)) {
        group_idx(root_idx_b) = root_idx_a;
        group_size(root_idx_a) += group_size(root_idx_b);
        group_size(root_idx_b) = 0;
    } else {
        group_idx(root_idx_a) = root_idx_b;
        group_size(root_idx_b) += group_size(root_idx_a);
        group_size(root_idx_a) = 0;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_single_group(amrex::Array1D<int, 1, NumSpec>& group_idx) {

    // This function checks whether all isotopes are either in the LIG group
    // or in another single group.

    const auto LIG_root_idx = get_root_index(NSE_INDEX::He4_index, group_idx);

    // Consider NSE when all nuclei belongs to the same group

    for (int n = 1; n <= NumSpec; ++n) {
        if (get_root_index(n, group_idx) != LIG_root_idx) {
            return false;
        }
    }
    return true;
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fill_reaction_timescale(amrex::Array1D<T, 1, NSE_INDEX::NumNSERatePairs>& reaction_timescales,
                             const int current_rate_pair, const amrex::Real rho,
                             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                             const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
                             const amrex::Real t_s) {
    // This function fills in the reaction timescale of the reaction for the RatePairs
    // where the forward and reverse are close, i.e. in equilibrium.

    // Calculate the forward and reverse rates of the current rate pair.

    constexpr amrex::Real inv6 = 1.0_rt / 6.0_rt;
    const auto r1_idx      = NSE_INDEX::rate_pair_data(current_rate_pair, 1);
    const auto r2_idx      = NSE_INDEX::rate_pair_data(current_rate_pair, 2);
    const auto r3_idx      = NSE_INDEX::rate_pair_data(current_rate_pair, 3);
    const auto p1_idx      = NSE_INDEX::rate_pair_data(current_rate_pair, 4);
    const auto p2_idx      = NSE_INDEX::rate_pair_data(current_rate_pair, 5);
    const auto p3_idx      = NSE_INDEX::rate_pair_data(current_rate_pair, 6);
    const auto nonNHA1_idx = NSE_INDEX::rate_pair_data(current_rate_pair, 7);
    const auto nonNHA2_idx = NSE_INDEX::rate_pair_data(current_rate_pair, 8);
    const auto forward_idx = NSE_INDEX::rate_pair_data(current_rate_pair, 9);
    const auto reverse_idx = NSE_INDEX::rate_pair_data(current_rate_pair, 10);

    amrex::Real b_f = screened_rates(forward_idx) * Y(r1_idx);
    amrex::Real b_r = screened_rates(reverse_idx) * Y(p1_idx);

    b_f *= (r2_idx != -1) ? (Y(r2_idx) * rho) : 1.0_rt;
    b_f *= (r3_idx != -1) ? (Y(r3_idx) * rho) : 1.0_rt;

    b_r *= (p2_idx != -1) ? (Y(p2_idx) * rho) : 1.0_rt;
    b_r *= (p3_idx != -1) ? (Y(p3_idx) * rho) : 1.0_rt;

    // Four possible case for multiplicity after sorting: ABC, AAB, ABB, AAA

    if (r1_idx == r2_idx && r2_idx == r3_idx) {
        // Case AAA
        b_f *= inv6;
    } else if (r1_idx == r2_idx || (r3_idx != -1 && r2_idx == r3_idx)) {
        // Case AAB or ABB. Need to make sure r3_idx is not -1
        b_f *= 0.5_rt;
    }

    if (p1_idx == p2_idx && p2_idx == p3_idx) {
        // Case AAA
        b_r *= inv6;
    } else if (p1_idx == p2_idx || (p3_idx != -1 && p2_idx == p3_idx)) {
        // Case AAB or ABB. Need to make sure p3_idx is not -1
        b_r *= 0.5_rt;
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

    amrex::Real b_min = amrex::min(b_f, b_r);

    // If one of the rate is zero, abort early.
    // We already set the reaction_timescales to be MAX_TIMESCALE
    // b_min should not be zero past this point.

    if (b_min == 0.0_rt) {
        return;
    }

    // Pre-filtered RatePair will at least have 1 nuclei not in (n, H1, He4)

    amrex::Real Y_min = Y(nonNHA1_idx);
    if (nonNHA2_idx != -1) {
        Y_min = amrex::min(Y_min, Y(nonNHA2_idx));
    }
    amrex::Real t_i = Y_min / b_min;

    //
    // Condition for checking if forward and reverse rates are in equilibrium
    // Also whether the reaction time scale is smaller than sound crossing time
    // See Equation 17 in Kushnir Paper
    //

    if ((std::abs(b_f - b_r) < 0.5_rt * ase_tol * (b_f + b_r)) &&
        (t_i < ase_tol * t_s)) {
        reaction_timescales(current_rate_pair) = t_i;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fill_merge_indices(int& merge_idx1, int& merge_idx2,
                        const int current_rate_pair,
                        amrex::Array1D<int, 1, NumSpec>& group_idx) {
    // Determines the merge indices for the current rate pair.
    // We consider three cases for each rate pair:
    //
    // 1. There are exactly three distinct groups: LIG and TWO other nuclei
    //    from TWO different groups. In this case, merge those two nuclei
    //    into the same group.
    //
    // 2. There are exactly TWO distinct groups: LIG and ONE or TWO other nuclei
    //    but belong to a SINGLE groups. In this case, merge those two groups
    //    into the same group.
    //
    // 3. There is exactly ONE distinct group: LIG -- all nuclei merged.
    //    Abort in this case.

    const auto nonNHA1_idx = NSE_INDEX::rate_pair_data(current_rate_pair, 7);
    const auto nonNHA2_idx = NSE_INDEX::rate_pair_data(current_rate_pair, 8);

    // We know that each RatePair will only have 1 or 2 non NHA nuclei.
    // First check if these nuclei are already merged into LIG.

    const auto LIG_root_idx = get_root_index(NSE_INDEX::He4_index, group_idx);
    const auto nonNHA1_root_idx = get_root_index(nonNHA1_idx, group_idx);

    // Use the first non NHA nuclei as the candidate for merging
    // Currently it has two possibilities: its in LIG or non in LIG
    // 1. First non NHA nuclei already in LIG
    // 2. First non NHA nuclei not in LIG, then set 2nd merge index to be LIG

    merge_idx1 = nonNHA1_root_idx;
    if (nonNHA1_root_idx != LIG_root_idx) {
        merge_idx2 = LIG_root_idx;
    }

    // If there is a second non NHA nuclei,
    // check if it is in the same group as the first non NHA nuclei.
    // If not, use this as the second merge index.

    if (nonNHA2_idx != -1) {
        const auto nonNHA2_root_idx = get_root_index(nonNHA2_idx, group_idx);
        if (nonNHA2_root_idx != nonNHA1_root_idx) {
            merge_idx2 = nonNHA2_root_idx;
        }
    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_grouping(amrex::Array1D<int, 1, NumSpec>& group_idx, const amrex::Real rho,
                  const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                  const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates,
                  const amrex::Real t_s) {
    // Groups all the nuclei using group_idx which contains the node number

    // Create a 1D array to keep track of the group size.
    // We want to merge smaller group to larger group during nse_union.

    amrex::Array1D<int, 1, NumSpec> group_size;

    // Initialize group_idx with 1-based indexing.
    // Also, initially all nuclei form their own group -- size = 1.

    for (int n = 1; n <= NumSpec; ++n) {
        group_idx(n) = n;
        group_size(n) = 1;
    }

    // let n,p,a form the same group (LIG) initially
    // Since we require He4 to be in the net, merge n and H1 with He4

    if constexpr (NSE_INDEX::N_index != -1) {
        nse_union(NSE_INDEX::N_index, NSE_INDEX::He4_index, group_idx, group_size);
    }
    if constexpr (NSE_INDEX::H1_index != -1) {
        nse_union(NSE_INDEX::H1_index, NSE_INDEX::He4_index, group_idx, group_size);
    }

    // Create a 1D array to hold reaction_timescales and RatePair indices
    // If forward and reverse rate are close, then store the reaction timescale,
    // otherwise set timescale to MAX_TIMESCALE -- slowest.
    // RatePair indices will hold index to the corresponding RatePair.

    constexpr amrex::Real MAX_TIMESCALE = std::numeric_limits<amrex::Real>::max();
    amrex::Array1D<int, 1, NSE_INDEX::NumNSERatePairs> rate_pair_idxs;
    amrex::Array1D<amrex::Real, 1, NSE_INDEX::NumNSERatePairs> reaction_timescales;

    for (int n = 1; n <= NSE_INDEX::NumNSERatePairs; ++n) {
        rate_pair_idxs(n) = n;
        reaction_timescales(n) = MAX_TIMESCALE;
        fill_reaction_timescale(reaction_timescales, n, rho, Y,
                                screened_rates, t_s);
    }

    // Sort rate_pair_idx using reaction_timescales.
    // The sorted rate_pair_idx should correspond to RatePairs
    // from smallest (fastest) to largest (slowest) timescale

    quickSort_Array1D(rate_pair_idxs, reaction_timescales);

    // After rate_pair_idx are sorted based on reaction timescales.
    // Do grouping based on timescale -- starting from the fastest timescale.

    for (int n = 1; n <= NSE_INDEX::NumNSERatePairs; ++n) {

        auto rate_pair_idx = rate_pair_idxs(n);

        // Check if the timescale is at max,
        // meaning reached the end of the valid RatePair -- break

        if (reaction_timescales(n) == MAX_TIMESCALE) {
            break;
        }

        // Fill in the nuclei indices for merging

        int merge_idx1 = -1;
        int merge_idx2 = -1;
        fill_merge_indices(merge_idx1, merge_idx2, rate_pair_idx, group_idx);

        // Check if merge indices are -1, which means nothing to merge

        if (merge_idx1 == -1 || merge_idx2 == -1) {
            continue;
        }

        // Merge nuclei into the same group

        nse_union(merge_idx1, merge_idx2, group_idx, group_size);
    }

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_nse(burn_t& current_state, bool skip_molar_check=false) {

    // This function returns the boolean that tells whether we're in nse or not
    // Note that it only works with pynucastro network for now.

    current_state.nse = false;

    amrex::Real T_in = current_state.T_fixed > 0.0_rt ? current_state.T_fixed : current_state.T;

    // If temperature is below T_min_nse
    // Or if we enable direct by a simple temperature threshold,
    // so T_nse_net > 0.0 but T_in < T_nse_net, we abort early to avoid computing cost.

    if ((T_in < T_min_nse) || (T_nse_net > 0.0_rt && T_in < T_nse_net)) {
        return current_state.nse;
    }

    // Check if any of the input rhoX is negative
    // This indicates a burn failure, so we should return false

    for (int n = 0; n < NumSpec; ++n) {
        if (current_state.y[SFS+n] < 0.0_rt) {
            return current_state.nse;
        }
    }

    // Get the NSE state using (rho, T, Ye) as input
    // Then compare input molar fractions to the NSE molar fractions.

    const auto nse_state = get_actual_nse_state(nse_input_rty, current_state);

    // Convert to molar fractions

    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    amrex::Array1D<amrex::Real, 1, NumSpec> Y_nse;
    const amrex::Real rho_inv = 1.0_rt / current_state.rho;

    for (int n = 0; n < NumSpec; ++n) {
        Y(n+1) = current_state.y[SFS+n] * aion_inv[n] * rho_inv;
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

    // Evaluate the rates.
    // Note that only some rates are needed, so perhaps we can optimize here
    // Or can we cache it from RHS calls during the integrator?

    rate_t rate_eval;
    constexpr int do_T_derivatives = 0;
    evaluate_rates<do_T_derivatives, rate_t>(state, Y, rate_eval);

    // need eos_state for speed of sound

    eos_t eos_state;

    // Initialize t_s, sound crossing timescale for a single zone to be max

    constexpr amrex::Real MAX_TIMESCALE = std::numeric_limits<amrex::Real>::max();
    amrex::Real t_s = MAX_TIMESCALE;

    // If we care about checking the reaction timescale to be smaller than t_s,

    if (!nse_dx_independent) {
        burn_to_eos(state, eos_state);
        eos(eos_input_rt, eos_state);
        t_s = state.dx / eos_state.cs;
    }

    // Create 1D array to hold the group node index

    amrex::Array1D<int, 1, NumSpec> group_idx;

    //
    // In Kushnir paper, we also need to perform a check
    // on finding at least 1 "fast reaction cycle"
    // However, it only works with neutron in the network.
    // I've also skipped this step when working with non-neutron
    // network. And it all seem to work fine.
    //

    // Now do nse grouping

    nse_grouping(group_idx, state.rho, Y, rate_eval.screened_rates, t_s);

    // Check if we result in a single group after grouping

    current_state.nse = in_single_group(group_idx);
    return current_state.nse;
}

#endif
