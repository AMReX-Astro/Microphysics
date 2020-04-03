#ifndef _vode_type_H_
#define _vode_type_H_

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <network.H>

const int VODE_NEQS = NumSpec + 2;

typedef amrex::Array1D<int, 1, VODE_NEQS> IArray1D;
typedef amrex::Array1D<Real, 1, VODE_NEQS> RArray1D;
typedef amrex::Array2D<Real, 1, VODE_NEQS, 1, VODE_NEQS> RArray2D;

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
    amrex::Real CONP, CRATE, DRC;
    amrex::Real ETA, ETAMAX, H, HNEW, HSCAL, PRL1;
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

    // Jacobian
    RArray2D jac;

    // Saved Jacobian
    RArray2D jac_save;

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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void vode_to_burn (const dvode_t& vode_state, burn_t& state)
{
    // Copy the integration data to the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        state.xn[n-1] = vode_state.y(n);
    }

    state.T = vode_state.y(net_itemp);
    state.e = vode_state.y(net_ienuc);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burn_to_vode (const burn_t& state, dvode_t& vode_state)
{
    // Copy the integration data from the burn state.

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.y(n) = state.xn[n-1];
    }

    vode_state.y(net_itemp) = state.T;
    vode_state.y(net_ienuc) = state.e;
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void renormalize_species (dvode_t& vode_state)
{
    Real sum = 0.0_rt;

    for (int n = 1; n <= NumSpec; ++n) {
        sum += vode_state.y(n);
    }

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.y(n) /= sum;
    }
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void clean_state (dvode_t& vode_state)
{

    // Ensure that mass fractions always stay positive and less than or
    // equal to 1.

    for (int n = 1; n <= NumSpec; ++n) {
        vode_state.y(n) = amrex::max(amrex::min(vode_state.y(n), 1.0_rt), SMALL_X_SAFE);
    }

    // Renormalize the abundances as necessary.

    if (renormalize_abundances) {
        renormalize_species(vode_state);
    }

    // Ensure that the temperature always stays within reasonable limits.

    vode_state.y(net_itemp) = amrex::min(MAX_TEMP, amrex::max(vode_state.y(net_itemp), EOSData::mintemp));

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_thermodynamics (burn_t& state, dvode_t& vode_state)
{

    // Fill an EOS state using a combination of the burn data
    // (which has supplementary data such as density) and the
    // VODE state (which has the current xn, T, and e).

    eos_t eos_state;

    burn_to_eos(state, eos_state);

    for (int n = 1; n <= NumSpec; ++n) {
        eos_state.xn[n-1] = vode_state.y(n);
    }
    eos_state.T = vode_state.y(net_itemp);
    eos_state.e = vode_state.y(net_ienuc);

    // Evaluate the thermodynamics -- if desired. Note that we can
    // either do this as an unconditional update (call_eos_in_rhs)
    // or we optionally set a fraction dT_crit such that we don't
    // call the EOS if the last temperature we evaluated the EOS
    // at is relatively close to the current temperature. Otherwise
    // just do the composition calculations since that's needed to
    // construct dX/dt.

    if (call_eos_in_rhs) {

        eos(eos_input_rt, eos_state);

    }
    else if (std::abs(eos_state.T - state.T_old) > dT_crit * eos_state.T)
    {

        eos(eos_input_rt, eos_state);

        state.dcvdT = (eos_state.cv - state.cv_old) / (eos_state.T - state.T_old);
        state.dcpdT = (eos_state.cp - state.cp_old) / (eos_state.T - state.T_old);

        state.T_old = eos_state.T;
        state.cv_old = eos_state.cv;
        state.cp_old = eos_state.cp;

    }
    else {

        composition(eos_state);

    }

    eos_to_burn(eos_state, state);

}

#endif
