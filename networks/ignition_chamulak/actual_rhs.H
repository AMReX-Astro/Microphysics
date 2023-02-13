#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>
#include <microphysics_math.H>
#include <actual_network.H>
#include <burn_type.H>
#include <rate_type.H>
#include <screen.H>
#include <jacobian_utilities.H>
#include <fundamental_constants.H>

using namespace amrex;
using namespace ArrayUtil;
using namespace Rates;

void actual_rhs_init ();

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real get_ebin (const Real& density)
{
    // Chamulak et al. provide the q-value resulting from C12 burning,
    // given as 3 different values (corresponding to 3 different densities).
    // Here we do a simple quadratic fit to the 3 values provided (see
    // Chamulak et al., p. 164, column 2).

    // our convention is that the binding energies are negative.  We convert
    // from the MeV values that are traditionally written in astrophysics
    // papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
    // MeV values are per nucleus, so we divide by aion to make it per
    // nucleon and we multiple by Avogardo's # (6.0221415e23) to get the
    // value in erg/g
    Real rho9 = density / 1.0e9_rt;

    // q_eff is effective heat evolved per reaction (given in MeV)
    Real q_eff = 0.06e0_rt * rho9 * rho9 + 0.02e0_rt * rho9 + 8.83e0_rt;

    // convert from MeV to ergs / gram.  Here M12_chamulak is the
    // number of C12 nuclei destroyed in a single reaction and 12.0 is
    // the mass of a single C12 nuclei.  Also note that our convention
    // is that binding energies are negative.
    return -q_eff * C::MeV2eV * C::ev2erg * C::n_A / (network::M12_chamulak * 12.0e0_rt);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rates (const burn_t& state, Real& rate, Real& dratedt, Real& sc1212, Real& dsc1212dt)
{
    using namespace Species;

    Real temp = state.T;
    Real dens = state.rho;

    // convert mass fractions to molar fractions
    Array1D<Real, 1, NumSpec> y;
    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // call the screening routine
    plasma_state_t pstate;
    fill_plasma_state(pstate, temp, dens, y);

    int jscr = 0;

    Real dsc1212dd;
    screen(pstate, jscr,
           zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1],
           sc1212, dsc1212dt, dsc1212dd);

    // compute some often used temperature constants
    Real T9     = temp / 1.e9_rt;
    Real dT9dt  = 1.0_rt / 1.e9_rt;
    Real T9a    = T9 / (1.0_rt + 0.0396e0_rt * T9);
    Real dT9adt = (T9a / T9 - (T9a / (1.0_rt + 0.0396e0_rt * T9)) * 0.0396e0_rt) * dT9dt;

    // compute the CF88 rate
    Real scratch    = std::pow(T9a, 1.0_rt / 3.0_rt);
    Real dscratchdt = (1.0_rt / 3.0_rt) * std::pow(T9a, -2.0_rt / 3.0_rt) * dT9adt;

    Real a    = 4.27e26_rt * std::pow(T9a, 5.0_rt / 6.0_rt) * std::pow(T9, -1.5e0_rt);
    Real dadt = (5.0_rt / 6.0_rt) * (a / T9a) * dT9adt - 1.5e0_rt * (a / T9) * dT9dt;

    Real b    = std::exp(-84.165e0_rt / scratch - 2.12e-3_rt * T9 * T9 * T9);
    Real dbdt = (84.165e0_rt * dscratchdt / (scratch * scratch) - 3.0_rt * 2.12e-3_rt * T9 * T9 * dT9dt) * b;

    rate = a * b;
    dratedt = dadt * b + a * dbdt;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ener_gener_rate (const Real& dydt_C12, const Real& ebin_C12)
{
    using namespace Species;

    return dydt_C12 * aion[C12-1] * ebin_C12;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (const burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    using namespace Species;

    for (int i = 1; i <= NumSpec; ++i) {
        ydot(i) = 0.0_rt;
    }

    Real rate, dratedt, sc1212, dsc1212dt;
    get_rates(state, rate, dratedt, sc1212, dsc1212dt);

    Real dens = state.rho;

    // we come in with mass fractions -- convert to molar fractions
    Array1D<Real, 1, NumSpec> y;
    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }


    // The change in number density of C12 is
    // d(n12)/dt = - M12_chamulak * 1/2 (n12)**2 <sigma v>
    //
    // where <sigma v> is the average of the relative velocity times the
    // cross section for the reaction, and the factor accounting for the
    // total number of particle pairs has a 1/2 because we are
    // considering a reaction involving identical particles (see Clayton
    // p. 293).  Finally, the -M12_chamulak means that for each reaction,
    // we lose M12_chamulak C12 nuclei (for a single rate, C12+C12,
    // M12_chamulak would be 2.  In Chamulak et al. (2008), they say a
    // value of 2.93 captures the energetics from a larger network
    //
    // Switching over to mass fractions, using n = rho X N_A/A, where N_A is
    // Avogadro's number, and A is the mass number of the nucleon, we get
    //
    // d(X12)/dt = -M12_chamulak * 1/2 (X12)**2 rho N_A <sigma v> / A12
    //
    // The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

    Real xc12tmp = amrex::max(y(C12) * aion[C12-1], 0.e0_rt);
    ydot(C12) = -(1.0_rt / 12.0_rt) * 0.5_rt * network::M12_chamulak * dens * sc1212 * rate * xc12tmp * xc12tmp;
    ydot(O16) = 0.0_rt;
    ydot(Ash) = -ydot(C12);

    // Convert back to molar form

    for (int i = 1; i <= NumSpec; ++i) {
        ydot(i) *= aion_inv[i-1];
    }

    Real ebin_C12 = get_ebin(dens);

    ydot(net_ienuc) = ener_gener_rate(ydot(C12), ebin_C12);
}

template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac (burn_t& state, MatrixType& jac)
{
    using namespace Species;

    Real rate, dratedt, sc1212, dsc1212dt;
    get_rates(state, rate, dratedt, sc1212, dsc1212dt);

    // initialize
    jac.zero();

    Real temp = state.T;
    Real dens = state.rho;
    Real xc12tmp = amrex::max(state.xn[C12-1], 0.0_rt);

    // Carbon jacobian elements

    jac(C12, C12) = -2.0_rt * (1.0_rt / 12.0_rt) * network::M12_chamulak * 0.5_rt * dens * sc1212 * rate * xc12tmp;
    jac(Ash, C12) = -jac(C12, C12);

    // Add the temperature derivative df(y_i) / dT (we will convert to d / de later)

    jac(C12, net_ienuc) = -(1.0_rt / 12.0_rt) * network::M12_chamulak * 0.5_rt *
                           (dens * rate * xc12tmp * xc12tmp * dsc1212dt +
                            dens * sc1212 * xc12tmp * xc12tmp * dratedt);

    // Convert back to molar form

    for (int i = 1; i <= neqs; ++i) {
        for (int j = 1; j <= NumSpec; ++j) {
            jac(j,i) *= aion_inv[j-1];
        }
    }

    // Energy generation rate Jacobian elements with respect to species

    Real ebin_C12 = get_ebin(dens);

    for (int j = 1; j <= NumSpec; ++j) {
        jac(net_ienuc, j) = ener_gener_rate(jac(C12, j), ebin_C12);
    }

    // df(e) / dT

    jac(net_ienuc, net_ienuc) = ener_gener_rate(jac(C12, net_ienuc), ebin_C12);

    // Convert to d / de

    jac(C12, net_ienuc) = temperature_to_energy_jacobian(state, jac(C12, net_ienuc));
    jac(net_ienuc, net_ienuc) = temperature_to_energy_jacobian(state, jac(net_ienuc, net_ienuc));
}

// Compute and store the more expensive screening factors

AMREX_INLINE
void set_up_screening_factors ()
{
    using namespace Species;

    // note: it is critical that these are called in the exact order
    // that the screening calls are done in the RHS routine, since we
    // use that order in the screening

    int jscr = 0;
    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);
}

#endif
