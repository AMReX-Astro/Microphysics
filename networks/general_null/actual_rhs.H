#ifndef actual_rhs_H
#define actual_rhs_H

#include <extern_parameters.H>
#include <actual_network.H>
#include <jacobian_utilities.H>

using namespace network_rp;

void actual_rhs_init ();

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_rhs (const burn_t& state, amrex::Array1D<amrex::Real, 1, neqs>& ydot)
{
    for (int i=1; i<=neqs; ++i) {
      ydot(i) = 0.;
    }
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_jac (const burn_t& state, MatrixType& jac)
{
    // initialize
    jac.zero();
}

#endif
