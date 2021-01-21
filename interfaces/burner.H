#ifndef _burner_H_
#define _burner_H_

#include <burn_type.H>
#ifdef SIMPLIFIED_SDC
#include <integrator_sdc.H>
#endif
#ifdef STRANG
#include <integrator.H>
#endif

#include <ArrayUtilities.H>

#ifdef NSE_TABLE
#include <nse.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void burner (burn_t& state, Real dt)
{

#ifdef STRANG

    // Strang code path

#ifdef NSE_TABLE

    bool nse_check = in_nse(state);

    if (nse_check) {

        // call the table
        nse_burn(state, dt);

    } else {
        // burn as usual
        integrator(state, dt);

        // update the aux from the new X's
        set_nse_aux_from_X(state);

    }

#else
    integrator(state, dt);
#endif

#else
#ifdef SIMPLIFIED_SDC
    // simplified SDC code path

    // right now, we don't have NSE implemented, so just call the integrator

    integrator(state, dt);

#endif
#endif

}

#endif
