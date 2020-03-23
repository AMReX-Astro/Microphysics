#ifndef _actual_network_data_H_
#define _actual_network_data_H_

#include <network_properties.H>

using namespace Species;

namespace iso7
{
    extern AMREX_GPU_MANAGED amrex::Real bion[NSpec];
    extern AMREX_GPU_MANAGED amrex::Real mion[NSpec];
    extern AMREX_GPU_MANAGED amrex::Real wion[NSpec];
}

namespace Rates {
    enum NetworkRates {ircag = 0,
                       iroga,
                       ir3a,
                       irg3a,
                       ir1212,
                       ir1216,
                       ir1616,
                       iroag,
                       irnega,
                       irneag,
                       irmgga,
                       irmgag,
                       irsiga,
                       ircaag,
                       irtiga,
                       irsi2ni,
                       irni2si,
                       NRates};

    const int NGroups = 4;

    extern amrex::Vector<std::string> names;
};

void actual_network_data_init()
{
    // Set the binding energy of the element
    bion(He4)  = 28.29603e0_rt;
    bion(C12)  = 92.16294e0_rt;
    bion(O16)  = 127.62093e0_rt;
    bion(Ne20) = 160.64788e0_rt;
    bion(Mg24) = 198.25790e0_rt;
    bion(Si28) = 236.53790e0_rt;
    bion(Ni56) = 484.00300e0_rt;

    // Set the mass
    for (int i = 0; i < NSpec; ++i) {
        mion(i) = nion(i) * C::mn + zion(i) * (C::mp + C::me) - bion(i) * C::mev2gr;

        // Molar mass
        wion(i) = C::avo * mion(i);

        // Common approximation to molar mass
        wion(i) = aion(i);
    }

    // set the names of the reaction rates
    {
        using namespace Rates;
        names.resize(NRates);
        names(ircag)   = "rcag";
        names(iroga)   = "roga";
        names(ir3a)    = "r3a";
        names(irg3a)   = "rg3a";    // inverse rate
        names(ir1212)  = "r1212";
        names(ir1216)  = "r1216";
        names(ir1616)  = "r1616";
        names(iroag)   = "roag";
        names(irnega)  = "rnega";
        names(irneag)  = "rneag";
        names(irmgga)  = "rmgga";
        names(irmgag)  = "rmgag";
        names(irsiga)  = "rsiga";
        names(ircaag)  = "rcaag";
        names(irtiga)  = "rtiga";
        names(irsi2ni) = "rsi2ni";
        names(irni2si) = "rni2si";
    }
}

#endif
