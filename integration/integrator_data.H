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
const int INT_NEQS = NumSpec + 1;
#endif


// We will use this parameter to determine if a given species
// abundance is unreasonably small or large (each X must satisfy
// -failure_tolerance <= X <= 1.0 + failure_tolerance).
const Real species_failure_tolerance = 1.e-2_rt;

enum integrator_errors {
    IERR_SUCCESS = 1,
    IERR_BAD_INPUTS = -1,
    IERR_DT_UNDERFLOW = -2,
    IERR_SPRAD_CONVERGENCE = -3,
    IERR_TOO_MANY_STEPS = -4,
    IERR_TOO_MUCH_ACCURACY_REQUESTED = -5,
    IERR_CORRECTOR_CONVERGENCE = -6,
    IERR_LU_DECOMPOSITION_ERROR = -7,
    IERR_ENTERED_NSE = -100
};


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
    int_neqs = NumSpec + 1;
#endif

    return int_neqs;
}

typedef amrex::Array1D<short, 1, INT_NEQS> IArray1D;
typedef amrex::Array1D<Real, 1, INT_NEQS> RArray1D;
typedef ArrayUtil::MathArray2D<1, INT_NEQS, 1, INT_NEQS> RArray2D;

#endif
