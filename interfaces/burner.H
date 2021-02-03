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

#ifdef NSE_TABLE

    bool nse_check = in_nse(state);

    if (nse_check) {

        // call the table
#ifdef SIMPLIFIED_SDC
        sdc_nse_burn(state, dt);
#else
        nse_burn(state, dt);
#endif

    } else {
        // burn as usual
        integrator(state, dt);

#ifndef SIMPLIFIED_SDC
        // update the aux from the new X's this is not needed with the
        // SDC method, since we took care of that internally.
        set_nse_aux_from_X(state);
#endif

    }

#else
    // we did not compile with NSE
    integrator(state, dt);
#endif

}

#endif
