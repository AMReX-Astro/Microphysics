#ifndef BURNER_H
#define BURNER_H

#include <burn_type.H>
#include <integrator.H>

#include <ArrayUtilities.H>

#ifdef NSE_TABLE
#include <nse_table_check.H>
#endif
#ifdef NSE_NET
#include <nse_check.H>
#endif
#ifdef NSE
#ifdef STRANG
#include <nse_update_strang.H>
#endif
#ifdef SDC
#include <nse_update_sdc.H>
#endif
#endif

template <typename BurnT>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burner (BurnT& state, amrex::Real dt)
{

#ifdef NSE

    bool nse_check = in_nse(state);

    if (nse_check) {

        // call the table
#ifdef SDC
        sdc_nse_burn(state, dt);
#endif

    } else {
        // burn as usual
        integrator(state, dt);

#if !defined(SDC) && defined(AUX_THERMO)
        // update the aux from the new X's this is not needed with the
        // SDC method, since we took care of that internally.
        set_aux_comp_from_X(state);
#endif

        // if we failed our burn, it may have been because we tried to
        // enter NSE and required too many steps.  At this point, the
        // aux data, T, and e represent the result from the incomplete
        // burn.

        // replace dt with just the remaining integration time
        // left after the failure
        amrex::Real dt_remaining = amrex::max(dt - state.time, 0.0_rt);

        // we use a relaxed NSE criteria now to catch states that are
        // right on the edge of being in NSE
#ifdef NSE_TABLE
        if (in_nse(state, true) && !state.success && dt_remaining > 0.0) {
#else
        if (in_nse(state, nse_skip_molar) && !state.success && dt_remaining > 0.0) {
#endif

#ifndef AMREX_USE_GPU
            std::cout << "recovering burn failure in NSE, zone = (" << state.i << ", " << state.j << ", " << state.k << ")" << std::endl;
#endif

            // This will append to state.e the amount additional
            // energy released from adjusting to the new NSE state
#ifdef SDC
            sdc_nse_burn(state, dt_remaining);
#endif
        }
    }

#ifndef AMREX_USE_GPU
    if (! state.success) {
        std::cout << state << std::endl;
        std::cout << in_nse(state) << std::endl;
        amrex::Error("unsuccessful burn");
    }
#endif
#else
    // we did not compile with NSE
    integrator(state, dt);
#endif

}

#endif
