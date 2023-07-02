#ifndef FE_TYPE_H
#define FE_TYPE_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>
#include <integrator_data.H>
#include <network.H>

#ifdef NETWORK_SOLVER
#include <actual_matrix.H>
#endif

namespace microphysics {

namespace forward_euler {

// When checking the integration time to see if we're done,
// be careful with roundoff issues.

const Real timestep_safety_factor = 1.0e-12_rt;

template <int int_neqs>
struct fe_t {

    amrex::Real t;      // the starting time
    amrex::Real tout;   // the stopping time

    int n_step;
    int n_rhs;
    int n_jac;

    amrex::Real atol_spec;
    amrex::Real rtol_spec;

    amrex::Real atol_enuc;
    amrex::Real rtol_enuc;

    Array1D<Real, 1, int_neqs> y;
    ArrayUtil::MathArray2D<1, int_neqs, 1, int_neqs> jac;

    short jacobian_type;
};

} // namespace forward_euler

} // namespace microphysics

#endif
