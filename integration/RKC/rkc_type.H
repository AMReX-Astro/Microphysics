#ifndef RKC_TYPE_H
#define RKC_TYPE_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>
#include <network.H>

#include <integrator_data.H>

const int RKC_NEQS = INT_NEQS;

typedef amrex::Array1D<Real, 1, RKC_NEQS> RArray1D;
typedef amrex::Array1D<Real, 1, RKC_NEQS> RArray1D;

const amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();

// We will use this parameter to determine if a given species abundance
// is unreasonably small or large (each X must satisfy
// -failure_tolerance <= X <= 1.0 + failure_tolerance).
const Real rkc_failure_tolerance = 1.e-2_rt;


struct rkc_t {

    // The initial point of the integration.
   amrex::Real t;

   // The end of the interval of integration.
   amrex::Real tend;

    // Integration array
    RArray1D y;

    // Tolerances
    amrex::Real rtol_spec, atol_spec;
    amrex::Real rtol_enuc, atol_enuc;
#ifdef TRUE_SDC
    amrex::Real rtol_dens, atol_dens;
#endif

    // RKC had an info array with 4 elements

    // info(1) was to allow for intermediate solutions -- we do not
    // use this (we do the equivalent of info(1) = 1)

    // info(2) defined how to compute the spectral radius.  We
    // have RKC do it internally (this is equivalent to info(2) = 0

    // info(3) was to allow for a constant Jacobian.  We state that
    // the Jacobian is not constant (this is equivalent to info(3) = 0

    // info(4) told RKC whether ATOL was a scalar or an array.
    // We have our own tolerances

    // work array -- we break the original Fortran work array up into
    // 5 separate arrays, each of length RKC_NEQS.  We don't need the
    // 8 other pieces of data that were stored in the Fortran version.
    //
    // for reference, the Fortran version was:
    // work(*) contains information needed for interpolation,
    // continuation after a return, and working storage. Items
    // relevant here are:
    //
    // The last step taken, hlast, is work(1).
    // The current t is work(2).
    // The number of equations, neqn, is work(3).
    // The unit roundoff, uround, is work(4).
    // The square root of uround, sqrtu, is work(5).
    // The maximum step size, hmax, is work(6).
    // The base address for the solution is ptr1 = nint(work(7)).
    // The solution at t starts at ptr1.
    // The derivative of the solution at t starts at ptr2.
    // The solution at t-hlast starts at ptr3.
    // The derivative of the solution at t-hlast starts at ptr4.
    // The estimated dominant eigenvector starts at ptr4 + neqn.

    amrex::Real hmax;

    RArray1D yn;
    RArray1D fn;
    RArray1D vtemp1;
    RArray1D vtemp2;
    RArray1D sprad;

    // set IDID = 0 for initialization
    short idid;

    // the following a variables that were in the RKCDID common
    // block

    // number of function evaluations
    int nfe;

    // number of integration stesp
    int nsteps;

    // number of accepted steps
    int naccpt;

    // number of rejected steps
    int nrejct;

    // number of evaluations of F used to estimate the spectral radius
    int nfesig;

    // maximum number of stages used
    int maxm;

};

#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#endif
#ifdef STRANG
#include <integrator_type_strang.H>
#endif

#endif
