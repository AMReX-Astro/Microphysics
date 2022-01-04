#ifndef NONAKA_PLOT_H
#define NONAKA_PLOT_H

#include <iomanip>
#include <sstream>
#include <string>

#include <extern_parameters.H>
#include <burn_type.H>
#include <network.H>

const int FIELD_WIDTH=14;

// this supports outputting the integration state at each RHS
// evaluation for a single zone.  This is used for diagnostics
// and for making the "Nonaka plot" used for comparing Strang
// and SDC integration

AMREX_INLINE
void nonaka_init() {


    std::stringstream ss;
    ss << std::setw(FIELD_WIDTH) << "time" << " ";

#ifdef SIMPLIFIED_SDC
    ss << std::setw(FIELD_WIDTH) << "sdc-iter" << " ";
#endif

#ifdef STRANG
    ss << std::setw(FIELD_WIDTH) << "strang-half" << " ";
#endif

    ss << std::setw(FIELD_WIDTH) << "density" << " ";
    ss << std::setw(FIELD_WIDTH) << "temperature" << " ";

    for (int i = 0; i < NumSpec; ++i) {
      ss << std::setw(FIELD_WIDTH) << "X(" + short_spec_names_cxx[i] + ")" << " ";
    }

    for (int i = 0; i < NumSpec; ++i) {
      ss << std::setw(FIELD_WIDTH) << "dXdt(" + short_spec_names_cxx[i] + ")" << " ";
    }

    ss << std::setw(FIELD_WIDTH) << "de/dt";

    std::ofstream nf;
    nf.open(nonaka_file, std::ios::out);
    nf << ss.str() << std::endl;
    nf.close();
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void nonaka_rhs(const Real time, const burn_t& state, const YdotNetArray1D& ydot_react) {

    // state: the burn_t corresponding to the current state note:
    // state.time is relative to the start of the current burn call,
    // so we rely on state.reference_time to be set from the calling
    // code with the offset.

    // ydot_react -- this is assumed to be the ydot we get directly
    // from the network -- it is not converted to the form needed by
    // the integrator yet (e.g., for SDC), so we should index it using
    // the network's convention.

    if (state.level == nonaka_level &&
        state.i == nonaka_i && state.j == nonaka_j && state.k == nonaka_k) {

        // append current state to nonaka log
        // at the current simulation time

        std::ofstream nf;
        nf.open(nonaka_file, std::ios::app);


        Real simulation_time = time + state.reference_time;

        nf << std::setw(FIELD_WIDTH) << simulation_time << " ";

#ifdef SIMPLIFIED_SDC
        nf << std::setw(FIELD_WIDTH) << state.sdc_iter << " ";
#endif

#ifdef STRANG
        nf << std::setw(FIELD_WIDTH) << state.strang_half << " ";
#endif

        nf << std::setw(FIELD_WIDTH) << state.rho << " ";

        nf << std::setw(FIELD_WIDTH) << state.T << " ";


        // Mass fractions X

        for (int n = 0; n < NumSpec; ++n) {
            nf << std::setw(FIELD_WIDTH) << state.xn[n] << " ";
        }

        // Convert molar fraction rhs to mass fraction rhs dX/dt

        // ydot_react is 1-based
        for (int n = 1; n <= NumSpec; ++n) {
            nf << std::setw(FIELD_WIDTH) << ydot_react(n) * aion[n-1] << " ";
        }

        // energy

        nf << std::setw(FIELD_WIDTH) << ydot_react(net_ienuc) << std::endl;

        nf.close();
    }

}

#endif
