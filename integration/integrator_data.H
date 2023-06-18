#ifndef INTEGRATOR_DATA_H
#define INTEGRATOR_DATA_H

#include <burn_type.H>

// Define the size of the ODE system that VODE will integrate

#ifdef STRANG
const int INT_NEQS = NumSpec + 1;
#endif

#ifdef SIMPLIFIED_SDC
const int INT_NEQS = SVAR_EVOLVE;
#endif

#ifdef TRUE_SDC
const int INT_NEQS = NumSpec + 2;
#endif

template<typename BurnT>
constexpr int integrator_neqs ()
{
    int int_neqs = 0;

#ifdef STRANG
    if constexpr (has_xn<BurnT>::value) {
        int_neqs += NumSpec;
    }

    if constexpr (has_energy<BurnT>::value) {
        int_neqs += 1;
    }
#endif

#ifdef SIMPLIFIED_SDC
    int_neqs = SVAR_EVOLVE;
#endif

#ifdef TRUE_SDC
    int_neqs = NumSpec + 2;
#endif

    return int_neqs;
}

typedef amrex::Array1D<short, 1, INT_NEQS> IArray1D;
typedef amrex::Array1D<Real, 1, INT_NEQS> RArray1D;
typedef ArrayUtil::MathArray2D<1, INT_NEQS, 1, INT_NEQS> RArray2D;

#endif
