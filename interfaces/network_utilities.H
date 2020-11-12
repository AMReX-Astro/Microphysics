#ifndef network_utilities_H
#define network_utilities_H

#include <AMReX_REAL.H>

struct rhs_t {
    amrex::Real prefactor;
    int specindex1;
    int specindex2;
    int specindex3;
};

// Form a unique numerical identifier from a given
// (species, rate) combination.
template<int num_rates>
constexpr int rhs_rate(int species, int rate)
{
    return rate + (species - 1) * num_rates;
}

#endif
