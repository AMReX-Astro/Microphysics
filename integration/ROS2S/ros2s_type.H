#ifndef ROS2S_TYPE_H
#define ROS2S_TYPE_H

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <ArrayUtilities.H>
#include <network.H>
#include <integrator_data.H>

#ifdef STRANG
#include <integrator_type_strang.H>
#endif
#ifdef SDC
#include <integrator_type_sdc.H>
#endif

// When checking the integration time to see if we're done,
// be careful with roundoff issues.

const amrex::Real timestep_safety_factor = 1.0e-4_rt;

template <int int_neqs>
struct ros2s_t {

    amrex::Real t;      // the starting time
    amrex::Real tout;   // the stopping time
    amrex::Real dt;     // next internal timestep

    int n_step;
    int n_rhs;
    int n_jac;

    amrex::Real atol_spec;
    amrex::Real rtol_spec;

    amrex::Real atol_enuc;
    amrex::Real rtol_enuc;

    amrex::Array1D<amrex::Real, 1, int_neqs> y;
    amrex::Array1D<amrex::Real, 1, int_neqs> ynew;
    amrex::Array1D<amrex::Real, 1, int_neqs> ak1;
    amrex::Array1D<amrex::Real, 1, int_neqs> ak2;
    amrex::Array1D<amrex::Real, 1, int_neqs> work;

    ArrayUtil::MathArray2D<jac_t, 1, int_neqs, 1, int_neqs> jac;
    IArray1D pivot;

    short jacobian_type;
};

#endif
