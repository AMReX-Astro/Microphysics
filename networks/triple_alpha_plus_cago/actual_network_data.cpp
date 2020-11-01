#include <AMReX_Vector.H>
#include <actual_network.H>

namespace triple_alpha_plus_cago
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}

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

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);

        names[ir3a-1]   = "3agc";   //     3 He4 --> C12
        names[ircago-1] = "cago";   // C12 + He4 --> O16
    }
}
