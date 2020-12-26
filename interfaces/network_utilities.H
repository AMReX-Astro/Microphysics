#ifndef network_utilities_H
#define network_utilities_H

#include <AMReX_REAL.H>

struct rhs_t {
    amrex::Real prefactor;
    int specindex1;
    int specindex2;
    int specindex3;
};

struct jac_t {
    amrex::Real prefactor;
    int specindex1;
    int specindex2;
    int specindex3;

    amrex::Real rate_prefactor;
    int rate_specindex1;
    int rate_specindex2;
    int rate_specindex3;
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
