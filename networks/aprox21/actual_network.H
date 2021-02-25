#ifndef actual_network_H
#define actual_network_H

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>

#include <fundamental_constants.H>
#include <network_properties.H>

using namespace amrex;

void actual_network_init();

namespace C
{
    namespace Legacy
    {
        // These are the values of the constants used in the original aprox21
        constexpr amrex::Real m_n = 1.67492721184e-24_rt;
        constexpr amrex::Real m_p = 1.67262163783e-24_rt;
        constexpr amrex::Real m_e = 9.1093821545e-28_rt;

        constexpr amrex::Real eV2erg  = 1.60217648740e-12_rt;
        constexpr amrex::Real MeV2erg = eV2erg*1.0e6_rt;
        constexpr amrex::Real MeV2gr  = MeV2erg/(c_light*c_light);

        constexpr amrex::Real n_A = 6.0221417930e23_rt;

        // conversion factor for nuclear energy generation rate
        constexpr amrex::Real enuc_conv2 = -n_A * c_light * c_light;
    }
}

const std::string network_name = "aprox21";

namespace aprox21
{
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    extern AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
}

namespace Rates {
    enum NetworkRates {
                        ir3a   = 1,
                        irg3a  = 2,
                        ircag  = 3,
                        ir1212 = 4,
                        ir1216 = 5,
                        ir1616 = 6,
                        iroga  = 7,
                        iroag  = 8,
                        irnega = 9,
                        irneag = 10,
                        irmgga = 11,
                        irmgag = 12,
                        irsiga = 13,
                        irmgap = 14,
                        iralpa = 15,
                        iralpg = 16,
                        irsigp = 17,
                        irsiag = 18,
                        irsga  = 19,
                        irsiap = 20,
                        irppa  = 21,
                        irppg  = 22,
                        irsgp  = 23,
                        irsag  = 24,
                        irarga = 25,
                        irsap  = 26,
                        irclpa = 27,
                        irclpg = 28,
                        irargp = 29,
                        irarag = 30,
                        ircaga = 31,
                        irarap = 32,
                        irkpa  = 33,
                        irkpg  = 34,
                        ircagp = 35,
                        ircaag = 36,
                        irtiga = 37,
                        ircaap = 38,
                        irscpa = 39,
                        irscpg = 40,
                        irtigp = 41,
                        irtiag = 42,
                        ircrga = 43,
                        irtiap = 44,
                        irvpa  = 45,
                        irvpg  = 46,
                        ircrgp = 47,
                        ircrag = 48,
                        irfega = 49,
                        ircrap = 50,
                        irmnpa = 51,
                        irmnpg = 52,
                        irfegp = 53,
                        irfeag = 54,
                        irniga = 55,
                        irfeap = 56,
                        ircopa = 57,
                        ircopg = 58,
                        irnigp = 59,
                        ir52ng = 60,
                        ir53gn = 61,
                        ir53ng = 62,
                        ir54gn = 63,
                        irfepg = 64,
                        ircogp = 65,
                        irheng = 66,
                        irhegn = 67,
                        irhng  = 68,
                        irdgn  = 69,
                        irdpg  = 70,
                        irhegp = 71,
                        irpen   = 72,
                        irnep   = 73,
                        irn56ec = 74,
                        irpp    = 75,
                        ir33    = 76,
                        irhe3ag = 77,
                        ircpg  = 78,
                        irnpg  = 79,
                        ifa    = 80,
                        ifg    = 81,
                        iropg  = 82,
                        irnag  = 83,
                        ir54ng   = 84,
                        ir55gn   = 85,
                        ir55ng   = 86,
                        ir56gn   = 87,
                        irfe54ap = 88,
                        irco57pa = 89,
                        irfe56pg = 90,
                        irco57gp = 91,
                        irr1   = 92,
                        irs1   = 93,
                        irt1   = 94,
                        iru1   = 95,
                        irv1   = 96,
                        irw1   = 97,
                        irx1   = 98,
                        ir1f54 = 99,
                        ir2f54 = 100,
                        ir3f54 = 101,
                        ir4f54 = 102,
                        ir5f54 = 103,
                        ir6f54 = 104,
                        ir7f54 = 105,
                        ir8f54 = 106,
                        iralf1 = 107,
                        iralf2 = 108,
                        irfe56_aux1 = 109,
                        irfe56_aux2 = 110,
                        irfe56_aux3 = 111,
                        irfe56_aux4 = 112,
                        NumRates = irfe56_aux4
                      };

    const int NumGroups = 4;
};

#endif
