#ifndef network_utilities_H
#define network_utilities_H

#include <AMReX_REAL.H>

struct rhs_t {
    amrex::Real prefactor;
    int specindex1;
    int specindex2;
    int specindex3;
};

#endif
