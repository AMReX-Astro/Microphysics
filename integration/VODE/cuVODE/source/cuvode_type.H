#ifndef _cuvode_type_H_
#define _cuvode_type_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <cuvode_parameters.H>
#include <vode_rpar.H>

const amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();
const amrex::Real SRUR = std::sqrt(UROUND);
const amrex::Real CCMXJ = 0.2e0_rt;
const amrex::Real HMIN = 0.0_rt;
const amrex::Real HMXI = 0.0_rt;

// For the backward differentiation formula (BDF) integration
// the maximum order should be no greater than 5.
const int VODE_MAXORD = 5;
const int VODE_LMAX = VODE_MAXORD + 1;

// How many timesteps should pass before refreshing the Jacobian
const int max_steps_between_jacobian_evals = 50;

// Type dvode_t contains the integration solution and control variables
struct dvode_t
{
    amrex::Real HU;
    amrex::Real ACNRM, CONP, CRATE, DRC;
    amrex::Real ETA, ETAMAX, H, HNEW, HSCAL, PRL1;
    amrex::Real RC, RL1, tn;
    int NFE, NJE, NST;
    int ICF, IPUP, JCUR, JSTART, JSV, kflag;
    int L;
    int MXSTEP;
    int NEWH, NEWQ, NQ, NQNYH, NQWAIT, NSLJ;
    int NSLP;

    amrex::Array1D<Real, 1, 13> el;
    amrex::Array1D<Real, 1, 13> tau;
    amrex::Array1D<Real, 1, 5> tq;

    // Tolerances
    amrex::Array1D<Real, 1, VODE_NEQS> rtol;
    amrex::Array1D<Real, 1, VODE_NEQS> atol;

    // Real parameters
    amrex::Array1D<Real, 1, n_rpar_comps> rpar;

    // State flag
    int istate;

    // Local time and integration end time
    amrex::Real t, tout;

    // Integration array
    amrex::Array1D<Real, 1, VODE_NEQS> y;

    // Jacobian
    amrex::Array2D<Real, 1, VODE_NEQS, 1, VODE_NEQS> jac;

    // Saved Jacobian
    amrex::Array2D<Real, 1, VODE_NEQS, 1, VODE_NEQS> jac_save;

    amrex::Array2D<Real, 1, VODE_NEQS, 1, VODE_LMAX> yh;
    amrex::Array1D<Real, 1, VODE_NEQS> ewt;
    amrex::Array1D<Real, 1, VODE_NEQS> savf;
    amrex::Array1D<Real, 1, VODE_NEQS> acor;

    // Jacobian method
    int jacobian;
};

#endif

