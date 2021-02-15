#ifndef FORTRAN_TO_CXX_ACTUAL_RHS_H
#define FORTRAN_TO_CXX_ACTUAL_RHS_H

#ifndef NETWORK_HAS_CXX_IMPLEMENTATION

#include <network.H>
#include <burn_type.H>

extern "C" void
fortran_to_cxx_actual_rhs(Real rho, Real T, Real e, const Real* xn,
                          Real abar, Real zbar, Real y_e, Real eta,
#if NAUX_NET > 0
                          const Real* aux,
#endif
#ifndef SIMPLIFIED_SDC
                          Real cv, Real cp, int self_heat,
                          Real T_old, Real dcvdT, Real dcpdT,
#endif
                          Real* const ydot);

extern "C" void
fortran_to_cxx_actual_jac(Real rho, Real T, Real e, const Real* xn,
                          Real abar, Real zbar, Real y_e, Real eta,
#if NAUX_NET > 0
                          const Real* aux,
#endif
#ifndef SIMPLIFIED_SDC
                          Real cv, Real cp, int self_heat,
                          Real T_old, Real dcvdT, Real dcpdT,
#endif
                          Real* const jac);

void
AMREX_FORCE_INLINE
actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
#ifndef SIMPLIFIED_SDC
    int self_heat = state.self_heat ? 1 : 0;
#endif

    fortran_to_cxx_actual_rhs(state.rho, state.T, state.e, state.xn,
                              state.abar, state.zbar, state.y_e, state.eta,
#if NAUX_NET > 0
                              state.aux,
#endif
#ifndef SIMPLIFIED_SDC
                              state.cv, state.cp, self_heat,
                              state.T_old, state.dcvdT, state.dcpdT,
#endif
                              &ydot(1));
}

template<class MatrixType>
void
AMREX_FORCE_INLINE
actual_jac(burn_t& state, MatrixType& jac)
{
#ifndef SIMPLIFIED_SDC
    int self_heat = state.self_heat ? 1 : 0;
#endif

    fortran_to_cxx_actual_jac(state.rho, state.T, state.e, state.xn,
                              state.abar, state.zbar, state.y_e, state.eta,
#if NAUX_NET > 0
                              state.aux,
#endif
#ifndef SIMPLIFIED_SDC
                              state.cv, state.cp, self_heat,
                              state.T_old, state.dcvdT, state.dcpdT,
#endif
                              &jac(1,1));
}

#endif

#endif
