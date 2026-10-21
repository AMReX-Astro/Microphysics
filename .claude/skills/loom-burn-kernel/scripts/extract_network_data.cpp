// Host extractor for gen_network_data.py: evaluates the stock constexpr data of
// an rhs.H-framework network (NEW_NETWORK_IMPLEMENTATION) and prints it, one
// record per line, doubles as C99 hex floats (%a, exact).
//
// Compiled by gen_network_data.py with the build's own -D/-I/-std flags (from
// `make -n` or an inventory flags.txt), from the unit test directory, plus
//   -I<out>  for extract_names.inc, written by the driver:
//            MPL_X_SPECIES(name) for every Species:: enumerator,
//            MPL_X_RATE(name) for every Rates:: enumerator except NumRates
//   -DMPL_X_HAVE_BION_MION=0/1  whether actual_network.H defines
//            network::bion<spec>() / network::mion<spec>()
// Nothing here runs on a device; the values are what the stock headers fold
// at compile time (rhs_data, locate_intermediate_rate_index,
// calculate_screen_factor) or define as constants.
#include <cstddef>
#include <cstdio>

#include <AMReX_Loop.H>

#include <network.H>
#include <burn_type.H>
#include <fundamental_constants.H>
#include <rhs.H>
#ifdef SCREENING
#include <screen_data.H>
#endif

#ifndef NEW_NETWORK_IMPLEMENTATION
#error "extract_network_data: this network does not use networks/rhs.H (NEW_NETWORK_IMPLEMENTATION undefined)"
#endif

#ifndef MPL_X_HAVE_BION_MION
#define MPL_X_HAVE_BION_MION 0
#endif

namespace {

void pr_def(const char* name, int on) { std::printf("DEF %s %d\n", name, on); }
void pr_int(const char* name, long v) { std::printf("CFG %s %ld\n", name, v); }
void pr_real(const char* tag, const char* name, double v) { std::printf("%s %s %a\n", tag, name, v); }

template <int rate>
void print_rhs_data()
{
    constexpr RHS::rhs_t d = RHS::rhs_data(rate);
    std::printf("RD %d species %d %d %d %d %d %d\n", rate,
                d.species_A, d.species_B, d.species_C, d.species_D, d.species_E, d.species_F);
    std::printf("RD %d number %a %a %a %a %a %a\n", rate,
                d.number_A, d.number_B, d.number_C, d.number_D, d.number_E, d.number_F);
    std::printf("RD %d exponent %d %d %d %d %d %d\n", rate,
                d.exponent_A, d.exponent_B, d.exponent_C, d.exponent_D, d.exponent_E, d.exponent_F);
    std::printf("RD %d forward_branching_ratio %a\n", rate, d.forward_branching_ratio);
    std::printf("RD %d reverse_branching_ratio %a\n", rate, d.reverse_branching_ratio);
    std::printf("RD %d apply_identical_particle_factor %d\n", rate, d.apply_identical_particle_factor);
    std::printf("RD %d rate_can_be_tabulated %d\n", rate, d.rate_can_be_tabulated);
    std::printf("RD %d screen_forward_reaction %d\n", rate, d.screen_forward_reaction);
    std::printf("RD %d screen_reverse_reaction %d\n", rate, d.screen_reverse_reaction);
    std::printf("RD %d additional_reaction %d %d %d\n", rate,
                d.additional_reaction_1, d.additional_reaction_2, d.additional_reaction_3);
    constexpr int index = RHS::locate_intermediate_rate_index(rate);
    std::printf("INTER %d %d\n", rate, index);
    std::printf("DENSEXP %d %d %d\n", rate, RHS::density_exponent_forward<rate>(),
                RHS::density_exponent_reverse<rate>());
}

#ifdef SCREENING
// One screened pair, evaluated as networks/rhs.H apply_screening does
// (`constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2)`).
template <int rate, int k>
void print_pair(const scrn::screen_factors_t& s)
{
    std::printf("SCREEN %d %d %a %a %a %a", rate, k, s.z1, s.a1, s.z2, s.a2);
#if SCREEN_METHOD == SCREEN_METHOD_screen5
    std::printf(" zs53=%a z1_53=%a z2_53=%a zs13=%a zs13inv=%a zhat=%a zhat2=%a lzav=%a aznut=%a",
                s.zs53, s.z1_53, s.z2_53, s.zs13, s.zs13inv, s.zhat, s.zhat2, s.lzav, s.aznut);
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2007
    std::printf(" ztilde=%a ztilde3=%a", s.ztilde, s.ztilde3);
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2009
    std::printf(" zs52=%a z1_52=%a z2_52=%a zs53=%a z1_53=%a z2_53=%a aznut=%a ztilde=%a",
                s.zs52, s.z1_52, s.z2_52, s.zs53, s.z1_53, s.z2_53, s.aznut, s.ztilde);
#elif SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    std::printf(" zs53=%a z1_53=%a z2_53=%a zs13=%a zs13inv=%a aznut=%a",
                s.zs53, s.z1_53, s.z2_53, s.zs13, s.zs13inv, s.aznut);
#endif
    std::printf("\n");
}

// Same branch structure as networks/rhs.H apply_screening.
template <int rate>
void print_screening()
{
    constexpr RHS::rhs_t data = RHS::rhs_data(rate);

    if constexpr (data.screen_forward_reaction == 0 && data.screen_reverse_reaction == 0) {
        return;
    }

    if constexpr (data.exponent_A == 1 && data.exponent_B == 1 && data.exponent_C == 0) {
        constexpr amrex::Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr amrex::Real A1 = NetworkProperties::aion(data.species_A);
        constexpr amrex::Real Z2 = NetworkProperties::zion(data.species_B);
        constexpr amrex::Real A2 = NetworkProperties::aion(data.species_B);
        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);
        static_assert(scn_fac.z1 == Z1);
        print_pair<rate, 1>(scn_fac);
    }

    if constexpr (data.exponent_A == 2 && data.exponent_B == 0 && data.exponent_C == 0) {
        constexpr amrex::Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr amrex::Real A1 = NetworkProperties::aion(data.species_A);
        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z1, A1);
        static_assert(scn_fac.z1 == Z1);
        print_pair<rate, 1>(scn_fac);
    }

    if constexpr (data.exponent_A == 3 && data.exponent_B == 0 && data.exponent_C == 0) {
        constexpr amrex::Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr amrex::Real A1 = NetworkProperties::aion(data.species_A);
        constexpr auto scn_fac1 = scrn::calculate_screen_factor(Z1, A1, Z1, A1);
        static_assert(scn_fac1.z1 == Z1);
        print_pair<rate, 1>(scn_fac1);
        constexpr amrex::Real Z2 = 2.0_rt * Z1;
        constexpr amrex::Real A2 = 2.0_rt * A1;
        constexpr auto scn_fac2 = scrn::calculate_screen_factor(Z1, A1, Z2, A2);
        static_assert(scn_fac2.z1 == Z1);
        print_pair<rate, 2>(scn_fac2);
    }
}
#endif

} // namespace

int main()
{
    std::printf("NETWORK %s\n", std::string(network_name).c_str());

    // configuration defines the kernel layout depends on
#ifdef SCREENING
    pr_def("SCREENING", 1);
#else
    pr_def("SCREENING", 0);
#endif
#ifdef SCREEN_METHOD
#if SCREEN_METHOD == SCREEN_METHOD_screen5
    std::printf("SCREEN_METHOD screen5\n");
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2007
    std::printf("SCREEN_METHOD chugunov2007\n");
#elif SCREEN_METHOD == SCREEN_METHOD_chugunov2009
    std::printf("SCREEN_METHOD chugunov2009\n");
#elif SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    std::printf("SCREEN_METHOD chabrier1998\n");
#elif SCREEN_METHOD == SCREEN_METHOD_null
    std::printf("SCREEN_METHOD null\n");
#else
    std::printf("SCREEN_METHOD unknown\n");
#endif
#endif
#ifdef NEUTRINOS
    pr_def("NEUTRINOS", 1);
#else
    pr_def("NEUTRINOS", 0);
#endif
#ifdef NSE
    pr_def("NSE", 1);
#else
    pr_def("NSE", 0);
#endif
#ifdef NSE_NET
    pr_def("NSE_NET", 1);
#else
    pr_def("NSE_NET", 0);
#endif
#ifdef NSE_TABLE
    pr_def("NSE_TABLE", 1);
#else
    pr_def("NSE_TABLE", 0);
#endif
#ifdef NONAKA_PLOT
    pr_def("NONAKA_PLOT", 1);
#else
    pr_def("NONAKA_PLOT", 0);
#endif
#ifdef STRANG
    pr_def("STRANG", 1);
#else
    pr_def("STRANG", 0);
#endif
#ifdef SIMPLIFIED_SDC
    pr_def("SIMPLIFIED_SDC", 1);
#else
    pr_def("SIMPLIFIED_SDC", 0);
#endif
#ifdef TRUE_SDC
    pr_def("TRUE_SDC", 1);
#else
    pr_def("TRUE_SDC", 0);
#endif
#ifdef SINGLE_PRECISION_JACOBIAN
    pr_def("SINGLE_PRECISION_JACOBIAN", 1);
#else
    pr_def("SINGLE_PRECISION_JACOBIAN", 0);
#endif
    pr_def("HAVE_BION_MION", MPL_X_HAVE_BION_MION);
    pr_int("sizeof_Real", static_cast<long>(sizeof(amrex::Real)));

    // sizes and burn_t indices (network_properties.H, interfaces/burn_type.H)
    pr_int("NumSpec", NumSpec);
    pr_int("NumSpecExtra", NumSpecExtra);
    pr_int("NumSpecTotal", NumSpecTotal);
    pr_int("NumAux", NumAux);
    pr_int("NumRates", Rates::NumRates);
    pr_int("num_intermediate", RHS::num_intermediate_reactions());
    pr_int("neqs", neqs);
    pr_int("net_ienuc", net_ienuc);
    pr_int("SFS", SFS);
    pr_int("SEINT", SEINT);
    pr_int("SFX", SFX);
    pr_int("SRHO", SRHO);
    pr_int("SMX", SMX);
    pr_int("SMY", SMY);
    pr_int("SMZ", SMZ);
    pr_int("SEDEN", SEDEN);
    pr_int("SVAR", SVAR);
    pr_int("SVAR_EVOLVE", SVAR_EVOLVE);
    pr_int("sizeof_burn_t", static_cast<long>(sizeof(burn_t)));

    // enumerators (names from the driver, values from the headers)
#define MPL_X_SPECIES(name) std::printf("SPECIES %s %d\n", #name, static_cast<int>(Species::name));
#define MPL_X_RATE(name) std::printf("RATE %s %d\n", #name, static_cast<int>(Rates::name));
#include "extract_names.inc"
#undef MPL_X_SPECIES
#undef MPL_X_RATE

    // per-species data: NetworkProperties::aion/zion (constexpr functions) and
    // the aion[] / zion[] / aion_inv[] arrays, which can differ in principle
    for (int sp = 1; sp <= NumSpecTotal; ++sp) {
        std::printf("ION %d %a %a %a %a %a\n", sp, NetworkProperties::aion(sp), NetworkProperties::zion(sp),
                    aion[sp - 1], zion[sp - 1], aion_inv[sp - 1]);
    }
#if MPL_X_HAVE_BION_MION
    amrex::constexpr_for<1, NumSpec + 1>([&] (auto n) {
        constexpr int sp = n;
        constexpr amrex::Real b = network::bion<sp>();
        constexpr amrex::Real m = network::mion<sp>();
        std::printf("BIONMION %d %a %a\n", sp, b, m);
    });
#endif

    // constants the per-network glue (ener_gener_rate) commonly uses
    pr_real("CONST", "legacy_enuc_conv2", C::Legacy::enuc_conv2);
    pr_real("CONST", "legacy_n_A", C::Legacy::n_A);
    pr_real("CONST", "legacy_MeV2erg", C::Legacy::MeV2erg);
    pr_real("CONST", "legacy_MeV2gr", C::Legacy::MeV2gr);
    pr_real("CONST", "enuc_conv2", C::enuc_conv2);
    pr_real("CONST", "n_A", C::n_A);

    // rhs_data for every rate, intermediate indices, density exponents
    amrex::constexpr_for<1, Rates::NumRates + 1>([&] (auto n) {
        constexpr int rate = n;
        print_rhs_data<rate>();
    });

#ifdef SCREENING
    amrex::constexpr_for<1, Rates::NumRates + 1>([&] (auto n) {
        constexpr int rate = n;
        print_screening<rate>();
    });
#endif

    std::printf("END\n");
    return 0;
}
