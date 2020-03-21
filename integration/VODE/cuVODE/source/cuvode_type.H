#ifndef _cuvode_type_H_
#define _cuvode_type_H_

#include <AMReX_REAL.H>

const amrex::Real UROUND = std::numeric_limits<amrex::Real>::min(1.0_rt);
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
    amrex::Real ACNRM, CONP, CRATE, DRC, EL[13];
    amrex::Real ETA, ETAMAX, H, HNEW, HSCAL, PRL1;
    amrex::Real RC, RL1, TAU[13], TQ[5], TN;
    int NFE, NJE, NST;
    int ICF, IPUP, JCUR, JSTART, JSV, KFLAG;
    int L;
    int MXSTEP;
    int NEWH, NEWQ, NQ, NQNYH, NQWAIT, NSLJ;
    int NSLP;

    // Tolerances
    amrex::Real RTOL[VODE_NEQS];
    amrex::Real ATOL[VODE_NEQS];

    // Real parameters
    amrex::Real RPAR[n_rpar_comps];

    // State flag
    int ISTATE;

    // Local time and integration end time
    amrex::Real T, TOUT;

    // Integration array
    amrex::Real Y[VODE_NEQS];

    // Jacobian
    amrex::Real jac[VODE_NEQS*VODE_NEQS];

    // Saved Jacobian
    amrex::Real jac_save[VODE_NEQS*VODE_NEQS];

    amrex::Real yh[VODE_LMAX][VODE_NEQS];
    amrex::Real ewt[VODE_NEQS];
    amrex::Real savf[VODE_NEQS];
    amrex::Real acor[VODE_NEQS];

    // Jacobian method
    int jacobian;
};

#endif

