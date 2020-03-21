#ifndef _cuvode_type_H_
#define _cuvode_type_H_

const Real UROUND = std::numeric_limits<Real>::min(1.0_rt);
const Real SRUR = std::sqrt(UROUND);
const Real CCMXJ = 0.2e0_rt;
const Real HMIN = 0.0_rt;
const Real HMXI = 0.0_rt;

// For the backward differentiation formula (BDF) integration
// the maximum order should be no greater than 5.
const int VODE_MAXORD = 5;
const int VODE_LMAX = VODE_MAXORD + 1;

// How many timesteps should pass before refreshing the Jacobian
const int max_steps_between_jacobian_evals = 50;

// Type dvode_t contains the integration solution and control variables
struct dvode_t
{
    Real HU;
    Real ACNRM, CONP, CRATE, DRC, EL[13];
    Real ETA, ETAMAX, H, HNEW, HSCAL, PRL1;
    Real RC, RL1, TAU[13], TQ[5], TN;
    int NFE, NJE, NST;
    int ICF, IPUP, JCUR, JSTART, JSV, KFLAG;
    int L;
    int MXSTEP;
    int NEWH, NEWQ, NQ, NQNYH, NQWAIT, NSLJ;
    int NSLP;

    // Tolerances
    Real RTOL[VODE_NEQS];
    Real ATOL[VODE_NEQS];

    // Real parameters
    Real RPAR[n_rpar_comps];

    // State flag
    int ISTATE;

    // Local time and integration end time
    Real T, TOUT;

    // Integration array
    Real Y[VODE_NEQS];

    // Jacobian
    Real jac[VODE_NEQS*VODE_NEQS];

    // Saved Jacobian
    Real jac_save[VODE_NEQS*VODE_NEQS];

    Real yh[VODE_LMAX][VODE_NEQS];
    Real ewt[VODE_NEQS];
    Real savf[VODE_NEQS];
    Real acor[VODE_NEQS];

    // Jacobian method
    int jacobian;
};

#endif

