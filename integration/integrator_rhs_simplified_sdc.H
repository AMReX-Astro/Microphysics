#ifndef INTEGRATOR_RHS_SIMPLIFIED_SDC_H
#define INTEGRATOR_RHS_SIMPLIFIED_SDC_H

#include <network.H>
#include <burn_type.H>
#include <extern_parameters.H>

#include <integrator_type_simplified_sdc.H>
#include <actual_network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_rhs.H>
#endif
#ifdef NONAKA_PLOT
#include <nonaka_plot.H>
#endif

// The f_rhs routine provides the right-hand-side for the integration solver.
// This is a generic interface that calls the specific RHS routine in the
// network you're actually using.

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs(const Real time, burn_t& state, T& int_state, RArray1D& ydot,
         [[maybe_unused]] const bool in_jacobian=false)
{

    // update rho

    update_density_in_time(time, state, int_state);

    // ensure that the mass fractions are valid -- only int_state is
    // updated here

    clean_state(time, state, int_state);

    // convert to the burn_t -- this does an EOS call to get T
    // and populates the (burn_t) state

    int_to_burn(time, int_state, state);

    // make sure that the temperature is valid

    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int n = 1; n <= INT_NEQS; ++n) {
            ydot(n) = 0.0_rt;
        }

        return;

    }

    // call the specific network to get the RHS

    actual_rhs(state, ydot);

#ifdef NONAKA_PLOT
    if (! in_jacobian) {
        nonaka_rhs(time, state, ydot);
    }
#endif

    // apply fudge factor:

    if (react_boost > 0.0_rt) {
        // the correct size here is the size of the reactive state
        // RHS, not the SDC state RHS (neqs, not INT_NEQS)

        for (int n = 1; n <= neqs; ++n) {
            ydot(n) *= react_boost;
        }
    }

    // convert back to the form needed by the integrator -- this will
    // add the advective terms

    rhs_to_int(time, state, ydot);

}



// Analytical Jacobian

template<typename T, class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (const Real time, burn_t& state, T& int_state, MatrixType& pd)
{

    // update rho, rho*u, ... in the burn_t state

    update_density_in_time(time, state, int_state);

    // ensure that the mass fractions are valid -- only int_state is
    // updated here

    clean_state(time, state, int_state);

    // convert to the burn_t -- this does an EOS call to get T
    // and populates the (burn_t) state

    int_to_burn(time, int_state, state);


    if (state.T <= EOSData::mintemp || state.T >= MAX_TEMP) {

        for (int j = 1; j <= INT_NEQS; ++j) {
            for (int i = 1; i <= INT_NEQS; ++i) {
                pd(i,j) = 0.0_rt;
            }
        }

        return;

    }

    // Call the specific network routine to get the Jacobian.

    actual_jac(state, pd);

    // The Jacobian from the nets is in terms of dYdot/dY, but we want
    // it was dXdot/dX, so convert here.
    for (int n = 1; n <= NumSpec; n++) {
        for (int m = 1; m <= neqs; m++) {
            pd(n,m) = pd(n,m) * aion[n-1];
        }
    }

    for (int m = 1; m <= neqs; m++) {
        for (int n = 1; n <= NumSpec; n++) {
            pd(m,n) = pd(m,n) * aion_inv[n-1];
        }
    }

    // apply fudge factor:

    if (react_boost > 0.0_rt) {
        pd.mul(react_boost);
    }

    // The system we integrate has the form (rho X_k, rho e)

    // pd is now of the form:
    //
    //  SFS         / d(rho X1dot)/dX1  d(rho X1dit)/dX2 ... 1/cv d(rho X1dot)/dT \ //
    //              | d(rho X2dot)/dX1  d(rho X2dot)/dX2 ... 1/cv d(rho X2dot)/dT |
    //  SFS-1+nspec |   ...                                                       |
    //  SEINT       \ d(rho Edot)/dX1   d(rho Edot)/dX2  ... 1/cv d(rho Edot)/dT  /
    //
    //                   SFS                                         SEINT

    // now correct the species derivatives
    // this constructs dy/dX_k |_e = dy/dX_k |_T - e_{X_k} |_T dy/dT / c_v

    eos_re_t eos_state;
    eos_state.rho = state.rho;
    eos_state.T = state.T;
    eos_state.e = state.e;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state.xn[n];
    }
#ifdef AUX_THERMO
    // make the aux data consistent with the state X's
    set_aux_comp_from_X(eos_state);
#endif

    eos(eos_input_re, eos_state);

    eos_xderivs_t eos_xderivs = composition_derivatives(eos_state);

    for (int m = 1; m <= neqs; m++) {
        for (int n = 1; n <= NumSpec; n++) {
            pd(m, n) -= eos_xderivs.dedX[n-1] * pd(m, net_ienuc);
        }
    }

}

#endif
