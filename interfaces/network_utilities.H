#ifndef network_utilities_H
#define network_utilities_H

#include <AMReX_REAL.H>

struct rhs_t {
    int species_A;
    int species_B;
    int species_C;
    int species_D;
    int species_E;
    int species_F;

    int number_A;
    int number_B;
    int number_C;
    int number_D;
    int number_E;
    int number_F;

    int exponent_A;
    int exponent_B;
    int exponent_C;
    int exponent_D;
    int exponent_E;
    int exponent_F;

    amrex::Real forward_branching_ratio;
    amrex::Real reverse_branching_ratio;
};

// Form a unique numerical identifier from a given
// (species, rate) combination.
template<int num_rates>
constexpr int rhs_rate (int species, int rate)
{
    return rate + (species - 1) * num_rates;
}

// Same as above but for (species1, species2, rate).
template<int num_species, int num_rates>
constexpr int jac_rate (int species1, int species2, int rate)
{
    return rate + (species2 - 1) * num_rates + (species1 - 1) * num_species * num_rates;
}

#endif
