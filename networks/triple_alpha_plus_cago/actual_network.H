#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

const std::string network_name = "triple_alpha_plus_cago";

namespace triple_alpha_plus_cago
{
    inline AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
}

namespace Rates
{
    enum NetworkRates {
        ir3a   = 1,
        ircago = 2,
        NumRates = ircago
    };

    const int NumGroups = 2;
};

AMREX_INLINE
void actual_network_init ()
{
    using namespace Species;
    using namespace triple_alpha_plus_cago;

    // our convention is that binding energy is negative.  The following are
    // the binding energies per unit mass (erg / g) obtained by converting
    // the energies in MeV to erg then multiplying by (N_A / aion) where
    // N_A = 6.0221415e23 is Avogadro's number
    bion(He4)  = 28.29603_rt; // MeV / nucleus
    bion(C12)  = 92.16294_rt; // MeV / nucleus
    bion(O16)  = 127.62093_rt; // MeV / nucleus
    bion(Fe56) = 492.25389_rt; // MeV / nucleus
}

#endif
