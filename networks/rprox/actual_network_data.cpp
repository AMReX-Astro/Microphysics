#include <AMReX_Vector.H>
#include <actual_network.H>

namespace rprox
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> ebin;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}

void actual_network_init()
{
    using namespace Species;
    using namespace rprox;

    // Our convention is that binding energy is negative.  The
    // following are the binding energies in MeV.
    ebin(C12)  = 92.16279_rt;
    ebin(O14)  = 98.7325_rt;
    ebin(O15)  = 111.9569_rt;
    ebin(O16)  = 127.6207_rt;
    ebin(F17)  = 128.2211_rt;
    ebin(Mg22) = 168.5768_rt;
    ebin(S30)  = 243.6866_rt;
    ebin(Ni56) = 483.995_rt;
    ebin(He4)  = 28.29599_rt;
    ebin(H1)   = 0.0_rt;

    // convert to erg / g by multiplying by N_A / aion and converting to erg
    for (int i = 1; i <= NumSpec; ++i) {
        ebin(i) = -ebin(i) * C::Legacy::N_A * C::Legacy::MeV2erg / aion[i-1];
    }

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NumRates);

        names[irlambCNO-1] = "rlambdaCNO";
        names[irag15o-1]   = "rag15o";
        names[irr1-1]      = "rr1";
        names[irag16o-1]   = "rag16o";
        names[irpg16o-1]   = "rpg16o";
        names[irpg17f-1]   = "rpg17f";
        names[irgp17f-1]   = "rgp17f";
        names[irlambda2-1] = "rlambda2";
        names[irap14o-1]   = "rap14o";
        names[irs1-1]      = "rs1";
        names[irlambda1-1] = "rlambda1";
        names[ir3a-1]      = "r3a";
        names[irpg12c-1]   = "rpg12c";
        names[irwk14o-1]   = "wk14o";
        names[irwk17f-1]   = "wk17f";
        names[irwk15o-1]   = "wk15o";
        names[irLweak-1]   = "Lweak";
        names[irla2-1]     = "la2";
    }
}
