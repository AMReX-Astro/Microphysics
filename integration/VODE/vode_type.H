#ifndef _vode_type_H_
#define _vode_type_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>
#include <network.H>

#ifdef NETWORK_SOLVER
#include <actual_matrix.H>
#endif

#ifdef SIMPLIFIED_SDC
const int VODE_NEQS = SVAR_EVOLVE;
#else
const int VODE_NEQS = NumSpec + 2;
#endif

typedef amrex::Array1D<int, 1, VODE_NEQS> IArray1D;
typedef amrex::Array1D<Real, 1, VODE_NEQS> RArray1D;
typedef ArrayUtil::MathArray2D<1, VODE_NEQS, 1, VODE_NEQS> RArray2D;

const amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();
const amrex::Real CCMXJ = 0.2e0_rt;
const amrex::Real HMIN = 0.0_rt;

// We will use this parameter to determine if a given species abundance
// is unreasonably small or large (each X must satisfy
// -failure_tolerance <= X <= 1.0 + failure_tolerance).
const Real vode_failure_tolerance = 1.e-2_rt;

// For each species whose abundance is above a certain threshold, we do
// not allow its mass fraction to change by more than a certain amount
// in any integration step.
const Real vode_increase_change_factor = 2.0_rt;
const Real vode_decrease_change_factor = 0.5_rt;

// For the backward differentiation formula (BDF) integration
// the maximum order should be no greater than 5.
const int VODE_MAXORD = 5;
const int VODE_LMAX = VODE_MAXORD + 1;

// How many timesteps should pass before refreshing the Jacobian
const int max_steps_between_jacobian_evals = 50;

// Type dvode_t contains the integration solution and control variables
struct dvode_t
{
    amrex::Real CONP, CRATE, DRC;
    amrex::Real ETA, ETAMAX, H, HNEW, HSCAL, PRL1, HMXI;
    amrex::Real RC, RL1, tn;
    int NFE, NJE, NST;
    int ICF, IPUP, JCUR;
    int L;
    int NEWH, NEWQ, NQ, NQNYH, NQWAIT, NSLJ;
    int NSLP;

    amrex::Array1D<Real, 1, VODE_LMAX> el;
    amrex::Array1D<Real, 1, VODE_LMAX> tau;
    amrex::Array1D<Real, 1, 5> tq;

    // Tolerances
    RArray1D rtol, atol;

    // Local time and integration end time
    amrex::Real t, tout;

    // Integration array
    RArray1D y;

#ifdef NETWORK_SOLVER

    // Jacobian
    SparseMatrix jac;

#ifndef AMREX_USE_GPU
    // Saved Jacobian
    SparseMatrix jac_save;
#endif

#else

    // Jacobian
    RArray2D jac;

#ifndef AMREX_USE_GPU
    // Saved Jacobian
    RArray2D jac_save;
#endif

#endif

    amrex::Array2D<Real, 1, VODE_NEQS, 1, VODE_LMAX> yh;

    RArray1D ewt, savf, acor;
};

#ifndef AMREX_USE_CUDA
AMREX_FORCE_INLINE
void print_state(dvode_t& dvode_state)
{

    std::cout << "CONP = " << dvode_state.CONP << std::endl;
    std::cout << "CRATE = " << dvode_state.CRATE << std::endl;
    std::cout << "DRC = " << dvode_state.DRC << std::endl;
    for (int n = 1; n <= VODE_LMAX; ++n) {
        std::cout << "el(" << n << ") = " << dvode_state.el(n) << std::endl;
    }
    std::cout << "el(1) = " << dvode_state.el(1) << std::endl;
    std::cout << "el(2) = " << dvode_state.el(2) << std::endl;
    std::cout << "el(3) = " << dvode_state.el(3) << std::endl;
    std::cout << "el(4) = " << dvode_state.el(4) << std::endl;
    std::cout << "el(5) = " << dvode_state.el(5) << std::endl;
    std::cout << "el(6) = " << dvode_state.el(6) << std::endl;
    std::cout << "ETA = " << dvode_state.ETA << std::endl;
    std::cout << "ETAMAX = " << dvode_state.ETAMAX << std::endl;
    std::cout << "H = " << dvode_state.H << std::endl;
    std::cout << "HNEW = " << dvode_state.HNEW << std::endl;
    std::cout << "HSCAL = " << dvode_state.HSCAL << std::endl;
    std::cout << "PRL1 = " << dvode_state.PRL1 << std::endl;
    std::cout << "RC = " << dvode_state.RC << std::endl;
    std::cout << "RL1 = " << dvode_state.RL1 << std::endl;
    for (int n = 1; n <= VODE_LMAX; ++n) {
        std::cout << "tau(" << n << ") = " << dvode_state.tau(n) << std::endl;
    }
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
    std::cout << "L = " << dvode_state.L << std::endl;
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

#ifdef SIMPLIFIED_SDC
#include <vode_type_simplified_sdc.H>
#endif
#ifdef STRANG
#include <vode_type_strang.H>
#endif

#endif
