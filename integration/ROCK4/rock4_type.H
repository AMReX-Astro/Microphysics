#ifndef ROCK4_TYPE_H
#define ROCK4_TYPE_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>
#include <network.H>

#include <integrator_data.H>

typedef amrex::Array1D<Real, 1, INT_NEQS> RArray1D;

const amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();

// We will use this parameter to determine if a given species abundance
// is unreasonably small or large (each X must satisfy
// -failure_tolerance <= X <= 1.0 + failure_tolerance).
const Real rock4_failure_tolerance = 1.e-2_rt;


struct rock4_t {

    // The initial point of the integration.
   amrex::Real t;

    // timestep
    amrex::Real h;

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
    // 8 separate arrays, each of length INT_NEQS.

    amrex::Real hmax;

    RArray1D yn;
    RArray1D fn;
    RArray1D fnt;
    RArray1D yjm1;
    RArray1D yjm2;
    RArray1D yjm3;
    RArray1D yjm4;
    RArray1D sprad;

    int mdeg;

    amrex::Array1D<int, 1, 2> mp;

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
    int nreject;

    // number of evaluations of F used to estimate the spectral radius
    int nfesig;

    // maximum number of stages used
    int max_stages;

    // maximum value of the spectral radius (nearest int)
    int max_rho;

    // minimum value of the spectral radius (nearest int)
    int min_rho;

};

#ifdef SIMPLIFIED_SDC
#include <integrator_type_simplified_sdc.H>
#endif
#ifdef STRANG
#include <integrator_type_strang.H>
#endif

#endif
