#ifndef actual_rhs_H
#define actual_rhs_H

#include <extern_parameters.H>
#include <actual_network.H>
#include <temperature_integration.H>
#include <rate_type.H>
#include <screen.H>

using namespace amrex;

void actual_rhs_init ();

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(burn_t& state, Array1D<rate_t, 1, Rates::NumGroups>& rr)
{
    using namespace Species;

    Real temp = state.T;
    Real dens = state.rho;

    // screening wants molar fractions
    Array1D<Real, 1, NumSpec> y;
    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // call the screening routine
    plasma_state_t pstate;
    fill_plasma_state(pstate, temp, dens, y);

    Real sc1212, dsc1212dt, dsc1212dd;
    
    int jscr = 0;
    screen5(pstate, jscr,
            zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1],
            sc1212, dsc1212dt, dsc1212dd);

    // compute some often used temperature constants
    Real T9     = temp / 1.e9_rt;
    Real dT9dt  = 1.0_rt / 1.e9_rt;
    Real T9a    = T9 / (1.0_rt + 0.0396_rt * T9);
    Real dT9adt = (T9a / T9 - (T9a / (1.0_rt + 0.0396_rt * T9)) * 0.0396_rt) * dT9dt;

    // compute the CF88 rate
    Real scratch    = std::pow(T9a, (1.0_rt / 3.0_rt));
    Real dscratchdt = (1.0_rt / 3.0_rt) * std::pow(T9a, (-2.0_rt / 3.0_rt)) * dT9adt;

    Real a       = 4.27e26_rt * std::pow(T9a, (5.0_rt / 6.0_rt)) * std::pow(T9, -1.5e0_rt);
    Real dadt    = (5.0_rt / 6.0_rt) * (a / T9a) * dT9adt - 1.5e0_rt * (a / T9) * dT9dt;

    Real b       = std::exp(-84.165e0_rt / scratch - 2.12e-3_rt * T9 * T9 * T9);
    Real dbdt    = (84.165e0_rt * dscratchdt / (scratch * scratch) - 3.0_rt * 2.12e-3_rt * T9 * T9 * dT9dt) * b;

    Real rate    = a * b;
    Real dratedt = dadt * b + a * dbdt;

    // Save the rates

    rr(1).rates(1)  = sc1212 * rate;
    rr(2).rates(1)  = dsc1212dt * rate + sc1212 * dratedt;
}


// Computes the instantaneous energy generation rate

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate (const Real& dydt, Real& enuc)
{
    using namespace Species;
    using namespace C::Legacy;
    using namespace ignition_simple;

    // This is basically e = m c**2

    // Note that since we don't explicitly evolve Mg24
    // in this network, we need to explicitly add its
    // contribution in this routine. We can factor out
    // the common factor of dydt(C12), we just need to
    // account for a factor of aion(C12) / aion(Mg24)
    // for the second term to make the expression work.

    enuc = dydt * (mion(C12) - mion(Mg24) / 2) * enuc_conv2;
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    using namespace Species;

    Array1D<rate_t, 1, Rates::NumGroups> rr;
    evaluate_rates(state, rr);

    // Now get the data from the state.

    Real temp = state.T;
    Real dens = state.rho;

    Real rate = rr(1).rates(1);

    // The change in number density of C12 is
    // d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
    //
    // where <sigma v> is the average of the relative velocity times the cross
    // section for the reaction, and the factor accounting for the total number
    // of particle pairs has a 1/2 because we are considering a reaction involving
    // identical particles (see Clayton p. 293).  Finally, the -2 means that for
    // each reaction, we lose 2 carbon nuclei.
    //
    // The corresponding Mg24 change is
    // d(n24)/dt = + 1/2 (n12)**2 <sigma v>
    //
    // note that no factor of 2 appears here, because we create only 1 Mg nuclei.
    //
    // Switching over to mass fractions, using n = rho X N_A/A, where N_A is
    // Avogadro's number, and A is the mass number of the nucleon, we get
    //
    // d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
    //
    // d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
    //
    // these are equal and opposite.
    //
    // The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

    Real xc12tmp = amrex::max(state.xn[C12-1], 0.0_rt);
    ydot(C12)  = -(1.0_rt / 12.0_rt) * dens * rate * xc12tmp * xc12tmp;
    ydot(O16)  = 0.0_rt;
    ydot(Mg24) = (1.0_rt / 12.0_rt) * dens * rate * xc12tmp * xc12tmp;

    // Convert back to molar form

    for (int i = 1; i <= NumSpec; ++i) {
        ydot(i) = ydot(i) * aion_inv[i-1];
    }

    ener_gener_rate(ydot(C12), ydot(net_ienuc));

#ifdef STRANG
    // Do the temperature equation explicitly here

    ydot(net_itemp) = temperature_rhs(state, ydot(net_ienuc));
#endif
}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac (burn_t& state, MatrixType& jac)
{
    using namespace Species;

    Array1D<rate_t, 1, Rates::NumGroups> rr;
    evaluate_rates(state, rr);

    // Get data from the state

    Real dens     = state.rho;

    Real rate     = rr(1).rates(1);
    Real dratedt  = rr(2).rates(1);
    Real xc12tmp  = amrex::max(state.xn[C12-1], 0.0_rt);

    // initialize
    jac.zero();

    // carbon jacobian elements
    jac(C12, C12) = -(1.0_rt / 6.0_rt) * dens * rate * xc12tmp;
    jac(Mg24, C12) = -jac(C12, C12);

    // add the temperature derivatives: df(y_i) / dT
    jac(C12, net_itemp)  = -(1.0_rt / 12.0_rt) * (dens * dratedt * xc12tmp * xc12tmp);

    // Convert back to molar form
    // Note that the factor of 1/A cancels in the (C12,C12) Jacobian element,
    // so this conversion is necessarily only for the temperature derivative.
    jac(C12, net_itemp) *= aion_inv[C12-1];

    // Energy generation rate Jacobian elements with respect to species
    ener_gener_rate(jac(C12, C12), jac(net_ienuc, C12));

    // Jacobian elements with respect to temperature
    ener_gener_rate(jac(C12, net_itemp), jac(net_ienuc, net_itemp));

    temperature_jac(state, jac);

}

#endif
