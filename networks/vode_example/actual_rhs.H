#ifndef actual_rhs_H
#define actual_rhs_H

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>

using namespace amrex;

AMREX_INLINE
void actual_rhs_init () {
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    using namespace Species;

    // Note: we are assuming that A = 1 for all species so molar and mass fractions are the same

    ydot(Y1) = -0.04_rt * state.xn[Y1-1] + 1.e4_rt * state.xn[Y2-1] * state.xn[Y3-1];
    ydot(Y2) = 0.04_rt * state.xn[Y1-1] - 1.e4_rt * state.xn[Y2-1] * state.xn[Y3-1] - 3.e7_rt * state.xn[Y2-1] * state.xn[Y2-1];
    ydot(Y3) = 3.e7_rt * state.xn[Y2-1] * state.xn[Y2-1];

    // this network has no energy generation

    ydot(net_ienuc) = 0.0_rt;

}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac (burn_t& state, MatrixType& jac)
{
    using namespace Species;

    // initialize
    jac.zero();

    // carbon jacobian elements
    jac(Y1, Y1) = -0.04_rt;
    jac(Y1, Y2) = 1.e4_rt * state.xn[Y3-1];
    jac(Y1, Y3) = 1.e4_rt * state.xn[Y2-1];

    jac(Y2, Y1) = 0.04_rt;
    jac(Y2, Y2) = -1.e4_rt * state.xn[Y3-1] - 6.e7_rt * state.xn[Y2-1];
    jac(Y2, Y3) = -1.e4_rt * state.xn[Y2-1];

    jac(Y3, Y1) = 0.0_rt;
    jac(Y3, Y2) = 6.e7_rt * state.xn[Y2-1];
    jac(Y3, Y3) = 0.0_rt;

    // all the derivatives involving energy are 0
}

#endif
