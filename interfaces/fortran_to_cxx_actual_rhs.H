#ifndef FORTRAN_TO_CXX_ACTUAL_RHS_H
#define FORTRAN_TO_CXX_ACTUAL_RHS_H

#ifndef NETWORK_HAS_CXX_IMPLEMENTATION

#include <network.H>
#include <burn_type.H>
#include <ArrayUtilities.H>
#include <jacobian_utilities.H>

extern "C" void
fortran_to_cxx_actual_rhs(Real rho, Real T, Real e, const Real* xn,
                          Real abar, Real zbar, Real y_e, Real eta,
#if NAUX_NET > 0
                          const Real* aux,
#endif
                          Real* const ydot);

extern "C" void
fortran_to_cxx_actual_jac(Real rho, Real T, Real e, const Real* xn,
                          Real abar, Real zbar, Real y_e, Real eta,
#if NAUX_NET > 0
                          const Real* aux,
#endif
                          Real* const jac);

void
AMREX_FORCE_INLINE
actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    // Create a temporary RHS which is one larger than ours (since Fortran
    // integrates both energy and temperature).

    Array1D<Real, 1, neqs+1> ydot_f;

    fortran_to_cxx_actual_rhs(state.rho, state.T, state.e, state.xn,
                              state.abar, state.zbar, state.y_e, state.eta,
#if NAUX_NET > 0
                              state.aux,
#endif
                              &ydot_f(1));

    // Copy the species components

    for (int n = 1; n <= NumSpec; ++n) {
        ydot(n) = ydot_f(n);
    }

    // Copy the energy component
    // Note that in Fortran, the temperature component is in what C++ would call
    // the energy slot (neqs), and the energy slot is in location neqs + 1.

    ydot(neqs) = ydot_f(net_ienuc + 1);
}

void
AMREX_FORCE_INLINE
actual_jac(burn_t& state, ArrayUtil::MathArray2D<1, neqs, 1, neqs>& jac)
{
    // Create a temporary Jacobian which is one larger than ours (since Fortran
    // integrates both energy and temperature).

    ArrayUtil::MathArray2D<1, neqs+1, 1, neqs+1> jac_f;

    fortran_to_cxx_actual_jac(state.rho, state.T, state.e, state.xn,
                              state.abar, state.zbar, state.y_e, state.eta,
#if NAUX_NET > 0
                              state.aux,
#endif
                              &jac_f(1,1));

    // Copy the df(species) / d(species) components

    for (int n = 1; n <= NumSpec; ++n) {
        for (int m = 1; m <= NumSpec; ++m) {
            jac(n, m) = jac_f(n, m);
        }
    }

    // Copy the df(enuc) / d(species) components
    // Note that in Fortran, the temperature component is in what C++ would call
    // the energy slot (neqs), and the energy slot is in location neqs + 1.

    for (int n = 1; n <= NumSpec; ++n) {
        jac(net_ienuc, n) = jac_f(net_ienuc + 1, n);
    }

    // Convert the df(species) / dT components

    for (int n = 1; n <= NumSpec; ++n) {
        jac(n, net_ienuc) = temperature_to_energy_jacobian(state, jac_f(n, net_ienuc));
    }

    // Convert the df(enuc) / dT component

    jac(net_ienuc, net_ienuc) = temperature_to_energy_jacobian(state, jac_f(net_ienuc+1, net_ienuc));
}

#endif

#endif
