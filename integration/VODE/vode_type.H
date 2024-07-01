#ifndef VODE_TYPE_H
#define VODE_TYPE_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <ArrayUtilities.H>
#include <network.H>

#include <integrator_data.H>

constexpr amrex::Real UROUND = std::numeric_limits<amrex::Real>::epsilon();

// CCMXJ  = Threshold on DRC for updating the Jacobian
constexpr amrex::Real CCMXJ = 0.2e0_rt;

constexpr amrex::Real HMIN = 0.0_rt;

// For each species whose abundance is above a certain threshold, we do
// not allow its mass fraction to change by more than a certain amount
// in any integration step.
constexpr amrex::Real vode_increase_change_factor = 4.0_rt;
constexpr amrex::Real vode_decrease_change_factor = 0.25_rt;

// For the backward differentiation formula (BDF) integration
// the maximum order should be no greater than 5.
constexpr int VODE_MAXORD = 5;
constexpr int VODE_LMAX = VODE_MAXORD + 1;

// How many timesteps should pass before refreshing the Jacobian
constexpr int max_steps_between_jacobian_evals = 50;

// Type dvode_t contains the integration solution and control variables
template<int int_neqs>
struct dvode_t
{
    // CONP   = The saved value of TQ(5)
    amrex::Real CONP;

    // CRATE  = Estimated corrector convergence rate constant
    amrex::Real CRATE;

    // Relative change in H*RL1 since last DVJAC call
    amrex::Real DRC;

    // ETA    = Saved tentative ratio of new to old H
    amrex::Real ETA;

    // ETAMAX = Saved maximum value of ETA to be allowed
    amrex::Real ETAMAX;

    // H      = The step size
    amrex::Real H;

    // HSCAL  = Stepsize in scaling of YH array
    amrex::Real HSCAL;

    // PRL1   = The saved value of RL1
    amrex::Real PRL1;

    // HMXI   = Inverse of the maximum absolute value of H to be used.
    //          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
    amrex::Real HMXI;

    // RC     = Ratio of current H*RL1 to value on last DVJAC call
    amrex::Real RC;

    // RL1    = The reciprocal of the coefficient EL(1)
    amrex::Real RL1;

    // TN     = The independent variable, updated on each step taken
    amrex::Real tn;

    // n_rhs    = The number of f evaluations for the problem so far
    int n_rhs;

    // n_jac    = The number of Jacobian evaluations so far
    int n_jac;

    // n_step    = The number of steps taken for the problem so far
    int n_step;

    // ICF    = Integer flag for convergence failure in DVNLSD:
    //            0 means no failures
    //            1 means convergence failure with out of date Jacobian
    //                   (recoverable error)
    //            2 means convergence failure with current Jacobian or
    //                   singular matrix (unrecoverable error)
    short ICF;

    // IPUP   = Saved flag to signal updating of Newton matrix
    short IPUP;

    // JCUR   = Output flag from DVJAC showing Jacobian status:
    //            JCUR = 0 means J is not current
    //            JCUR = 1 means J is current
    short JCUR;

    // L      = Integer variable, NQ + 1, current order plus one
    short L;

    // NEWH   = Saved integer to flag change of H
    short NEWH;

    // NEWQ   = The method order to be used on the next step
    short NEWQ;

    // NQ     = Integer variable, the current integration method order
    short NQ;

    // NQWAIT = A counter controlling the frequency of order changes.
    //          An order change is about to be considered if NQWAIT = 1.
    short NQWAIT;

    // NSLJ   = The number of steps taken as of the last Jacobian update
    int NSLJ;

    // NSLP   = Saved value of n_step as of last Newton matrix update
    int NSLP;

    // jacobian_type = the type of Jacobian to use (1 = analytic, 2 = numerical)
    short jacobian_type;

    // EL     = Real array of integration coefficients.  See DVSET
    amrex::Array1D<amrex::Real, 1, VODE_LMAX> el;

    // TAU    = Real vector of past NQ step sizes, length 13
    amrex::Array1D<amrex::Real, 1, VODE_LMAX> tau;

    // TQ     = A real vector of length 5 in which DVSET stores constants
    //          used for the convergence test, the error test, and the
    //          selection of H at a new order.
    amrex::Array1D<amrex::Real, 1, 5> tq;

    // Tolerances
    amrex::Real rtol_spec, atol_spec;
    amrex::Real rtol_enuc, atol_enuc;

    // Local time and integration end time
    amrex::Real t, tout;

    // Integration array
    amrex::Array1D<amrex::Real, 1, int_neqs> y;

    // Jacobian
    ArrayUtil::MathArray2D<1, int_neqs, 1, int_neqs> jac;

#ifdef ALLOW_JACOBIAN_CACHING
    // Saved Jacobian
    ArrayUtil::MathArray2D<1, int_neqs, 1, int_neqs> jac_save;
#endif

    // the Nordsieck history array
    amrex::Array2D<amrex::Real, 1, int_neqs, 1, VODE_LMAX> yh;

    amrex::Array1D<amrex::Real, 1, int_neqs> ewt, savf;

    amrex::Array1D<short, 1, int_neqs> pivot;

    // Array of size NEQ used for the accumulated corrections on each
    // step, scaled in the output to represent the estimated local
    // error in Y on the last step.  This is the vector e in the
    // description of the error control.  It is defined only on a
    // successful return from DVODE.
    amrex::Array1D<amrex::Real, 1, int_neqs> acor;
};

#ifndef AMREX_USE_GPU
template <int int_neqs>
AMREX_FORCE_INLINE
void print_state(dvode_t<int_neqs>& dvode_state)
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
    std::cout << "n_rhs = " << dvode_state.n_rhs << std::endl;
    std::cout << "n_jac = " << dvode_state.n_jac << std::endl;
    std::cout << "n_step = " << dvode_state.n_step << std::endl;
    std::cout << "ICF = " << dvode_state.ICF << std::endl;
    std::cout << "IPUP = " << dvode_state.IPUP << std::endl;
    std::cout << "JCUR = " << dvode_state.JCUR << std::endl;
    std::cout << "L = " << dvode_state.L << std::endl;
    std::cout << "NEWH = " << dvode_state.NEWH << std::endl;
    std::cout << "NEWQ = " << dvode_state.NEWQ << std::endl;
    std::cout << "NQ = " << dvode_state.NQ << std::endl;
    std::cout << "NQWAIT = " << dvode_state.NQWAIT << std::endl;
    std::cout << "NSLJ = " << dvode_state.NSLJ << std::endl;
    std::cout << "NSLP = " << dvode_state.NSLP << std::endl;

    for (int i = 1; i <= int_neqs; ++i) {
        std::cout << "y(" << i << ") = " << dvode_state.y(i) << std::endl;
    }

    for (int j = 1; j <= VODE_LMAX; ++j) {
        for (int i = 1; i <= int_neqs; ++i) {
            std::cout << "yh(" << i << "," << j << ") = " << dvode_state.yh(i,j) << std::endl;
        }
    }

    for (int i = 1; i <= int_neqs; ++i) {
        std::cout << "ewt(" << i << ") = " << dvode_state.ewt(i) << std::endl;
    }

    for (int i = 1; i <= int_neqs; ++i) {
        std::cout << "savf(" << i << ") = " << dvode_state.savf(i) << std::endl;
    }

    for (int i = 1; i <= int_neqs; ++i) {
        std::cout << "acor(" << i << ") = " << dvode_state.acor(i) << std::endl;
    }
}
#endif

#ifdef SDC
#include <integrator_type_sdc.H>
#endif
#ifdef STRANG
#include <integrator_type_strang.H>
#endif

#endif
