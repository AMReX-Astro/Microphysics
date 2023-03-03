#ifndef BE_TYPE_H
#define BE_TYPE_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>

#include <integrator_data.H>
#ifdef STRANG
#include <integrator_type_strang.H>
#endif
#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#endif
#include <network.H>

#ifdef NETWORK_SOLVER
#include <actual_matrix.H>
#endif

const int BE_NEQS = INT_NEQS;

const int BE_SUCCESS = 0;
const int BE_NONCONVERGENCE = -100;
const int BE_LU_DECOMPOSITION_ERROR = -101;
const int BE_TOO_MANY_STEPS = -102;

// When checking the integration time to see if we're done,
// be careful with roundoff issues.

const Real timestep_safety_factor = 1.0e-12_rt;

#define VODELIKE_ERROR 1

struct be_t {

    amrex::Real t;      // the starting time
    amrex::Real tout;   // the stopping time

    int n_step;
    int n_rhs;
    int n_jac;

    amrex::Real atol_spec;
    amrex::Real rtol_spec;

    amrex::Real atol_enuc;
    amrex::Real rtol_enuc;

    Array1D<Real, 1, BE_NEQS> y;
    JacNetArray2D jac;

    short jacobian_type;
};

#endif
