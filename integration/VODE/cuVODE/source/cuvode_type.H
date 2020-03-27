#ifndef _cuvode_type_H_
#define _cuvode_type_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <cuvode_parameters.H>
#include <vode_rpar.H>

const amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();
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

#ifndef AMREX_USE_CUDA
void print_state(dvode_t& dvode_state)
{

    std::cout << "HU = " << dvode_state.HU << std::endl;
    std::cout << "ACNRM = " << dvode_state.ACNRM << std::endl;
    std::cout << "CONP = " << dvode_state.CONP << std::endl;
    std::cout << "CRATE = " << dvode_state.CRATE << std::endl;
    std::cout << "DRC = " << dvode_state.DRC << std::endl;
    std::cout << "el(1) = " << dvode_state.el(1) << std::endl;
    std::cout << "el(2) = " << dvode_state.el(2) << std::endl;
    std::cout << "el(3) = " << dvode_state.el(3) << std::endl;
    std::cout << "el(4) = " << dvode_state.el(4) << std::endl;
    std::cout << "el(5) = " << dvode_state.el(5) << std::endl;
    std::cout << "el(6) = " << dvode_state.el(6) << std::endl;
    std::cout << "el(7) = " << dvode_state.el(7) << std::endl;
    std::cout << "el(8) = " << dvode_state.el(8) << std::endl;
    std::cout << "el(9) = " << dvode_state.el(9) << std::endl;
    std::cout << "el(10) = " << dvode_state.el(10) << std::endl;
    std::cout << "el(11) = " << dvode_state.el(11) << std::endl;
    std::cout << "el(12) = " << dvode_state.el(12) << std::endl;
    std::cout << "el(13) = " << dvode_state.el(13) << std::endl;
    std::cout << "ETA = " << dvode_state.ETA << std::endl;
    std::cout << "ETAMAX = " << dvode_state.ETAMAX << std::endl;
    std::cout << "H = " << dvode_state.H << std::endl;
    std::cout << "HNEW = " << dvode_state.HNEW << std::endl;
    std::cout << "HSCAL = " << dvode_state.HSCAL << std::endl;
    std::cout << "PRL1 = " << dvode_state.PRL1 << std::endl;
    std::cout << "RC = " << dvode_state.RC << std::endl;
    std::cout << "RL1 = " << dvode_state.RL1 << std::endl;
    std::cout << "tau(1) = " << dvode_state.tau(1) << std::endl;
    std::cout << "tau(2) = " << dvode_state.tau(2) << std::endl;
    std::cout << "tau(3) = " << dvode_state.tau(3) << std::endl;
    std::cout << "tau(4) = " << dvode_state.tau(4) << std::endl;
    std::cout << "tau(5) = " << dvode_state.tau(5) << std::endl;
    std::cout << "tau(6) = " << dvode_state.tau(6) << std::endl;
    std::cout << "tau(7) = " << dvode_state.tau(7) << std::endl;
    std::cout << "tau(8) = " << dvode_state.tau(8) << std::endl;
    std::cout << "tau(9) = " << dvode_state.tau(9) << std::endl;
    std::cout << "tau(10) = " << dvode_state.tau(10) << std::endl;
    std::cout << "tau(11) = " << dvode_state.tau(11) << std::endl;
    std::cout << "tau(12) = " << dvode_state.tau(12) << std::endl;
    std::cout << "tau(13) = " << dvode_state.tau(13) << std::endl;
    std::cout << "tq(1) = " << dvode_state.tq(1) << std::endl;
    std::cout << "tq(2) = " << dvode_state.tq(2) << std::endl;
    std::cout << "tq(3) = " << dvode_state.tq(3) << std::endl;
    std::cout << "tq(4) = " << dvode_state.tq(4) << std::endl;
    std::cout << "tq(5) = " << dvode_state.tq(5) << std::endl;
    std::cout << "tn = " << dvode_state.tn << std::endl;
    std::cout << "NFE = " << dvode_state.NFE << std::endl;
    std::cout << "NJE = " << dvode_state.NJE << std::endl;
    std::cout << "NST = " << dvode_state.NST << std::endl;
    std::cout << "ICF = " << dvode_state.ICF << std::endl;
    std::cout << "IPUP = " << dvode_state.IPUP << std::endl;
    std::cout << "JCUR = " << dvode_state.JCUR << std::endl;
    std::cout << "JSTART = " << dvode_state.JSTART << std::endl;
    std::cout << "JSV = " << dvode_state.JSV << std::endl;
    std::cout << "kflag = " << dvode_state.kflag << std::endl;
    std::cout << "L = " << dvode_state.L << std::endl;
    std::cout << "MXSTEP = " << dvode_state.MXSTEP << std::endl;
    std::cout << "NEWH = " << dvode_state.NEWH << std::endl;
    std::cout << "NEWQ = " << dvode_state.NEWQ << std::endl;
    std::cout << "NQ = " << dvode_state.NQ << std::endl;
    std::cout << "NQNYH = " << dvode_state.NQNYH << std::endl;
    std::cout << "NQWAIT = " << dvode_state.NQWAIT << std::endl;
    std::cout << "NSLJ = " << dvode_state.NSLJ << std::endl;
    std::cout << "NSLP = " << dvode_state.NSLP << std::endl;

    for (int i = 1; i <= VODE_NEQS; ++i) {
        std::cout << "y(" << i << ") = " << dvode_state.y(i) << std::endl;
    }

    for (int j = 1; j <= VODE_LMAX; ++j) {
        for (int i = 1; i <= VODE_NEQS; ++i) {
            std::cout << "yh(" << i << "," << j << ") = " << dvode_state.yh(i,j) << std::endl;
        }
    }

    for (int i = 1; i <= VODE_NEQS; ++i) {
        std::cout << "ewt(" << i << ") = " << dvode_state.ewt(i) << std::endl;
    }

    for (int i = 1; i <= VODE_NEQS; ++i) {
        std::cout << "savf(" << i << ") = " << dvode_state.savf(i) << std::endl;
    }

    for (int i = 1; i <= VODE_NEQS; ++i) {
        std::cout << "acor(" << i << ") = " << dvode_state.acor(i) << std::endl;
    }

}
#endif

#endif

