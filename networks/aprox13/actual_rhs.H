#ifndef _actual_rhs_H_
#define _actual_rhs_H_

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>
#include <microphysics_math.H>
#include <actual_network.H>
#include <burn_type.H>
#include <tfactors.H>
#include <screen.H>
#include <sneut5.H>
#include <aprox_rates.H>
#include <temperature_integration.H>

using namespace amrex;
using namespace ArrayUtil;

void actual_rhs_init();

namespace RateTable
{
    constexpr Real tab_tlo = 6.0e0_rt;
    constexpr Real tab_thi = 10.0e0_rt;
    constexpr int tab_per_decade = 500;
    constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
    constexpr Real tab_tstp = (tab_thi - tab_tlo) / static_cast<Real>(tab_imax - 1);

    extern AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> rattab;
    extern AMREX_GPU_MANAGED Array2D<Real, 1, Rates::NumRates, 1, nrattab> drattabdt;
    extern AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr int density_power(const int id)
    {

        using namespace Rates;

        if (id == ir3a) return 2;
        else if (id == irg3a) return 0;
        else if (id == ircag) return 1;
        else if (id == iroga) return 0;
        else if (id == ir1212) return 1;
        else if (id == ir1216) return 1;
        else if (id == ir1616) return 1;
        else if (id == iroag) return 1;
        else if (id == irnega) return 0;
        else if (id == irneag) return 1;
        else if (id == irmgga) return 0;
        else if (id == irmgag) return 1;
        else if (id == irsiga) return 0;
        else if (id == irmgap) return 1;
        else if (id == iralpa) return 1;
        else if (id == iralpg) return 1;
        else if (id == irsigp) return 0;
        else if (id == irsiag) return 1;
        else if (id == irsga) return 0;
        else if (id == irppa) return 1;
        else if (id == irsiap) return 1;
        else if (id == irppg) return 1;
        else if (id == irsgp) return 0;
        else if (id == irsag) return 1;
        else if (id == irarga) return 0;
        else if (id == irsap) return 1;
        else if (id == irclpa) return 1;
        else if (id == irclpg) return 1;
        else if (id == irargp) return 0;
        else if (id == irarag) return 1;
        else if (id == ircaga) return 0;
        else if (id == irarap) return 1;
        else if (id == irkpa) return 1;
        else if (id == irkpg) return 1;
        else if (id == ircagp) return 0;
        else if (id == ircaag) return 1;
        else if (id == irtiga) return 0;
        else if (id == ircaap) return 1;
        else if (id == irscpa) return 1;
        else if (id == irscpg) return 1;
        else if (id == irtigp) return 0;
        else if (id == irtiag) return 1;
        else if (id == ircrga) return 0;
        else if (id == irtiap) return 1;
        else if (id == irvpa) return 1;
        else if (id == irvpg) return 1;
        else if (id == ircrgp) return 0;
        else if (id == ircrag) return 1;
        else if (id == irfega) return 0;
        else if (id == ircrap) return 1;
        else if (id == irmnpa) return 1;
        else if (id == irmnpg) return 1;
        else if (id == irfegp) return 0;
        else if (id == irfeag) return 1;
        else if (id == irniga) return 0;
        else if (id == irfeap) return 1;
        else if (id == ircopa) return 1;
        else if (id == ircopg) return 1;
        else if (id == irnigp) return 0;
        else if (id == irr1) return -1;
        else if (id == irs1) return -1;
        else if (id == irt1) return -1;
        else if (id == iru1) return -1;
        else if (id == irv1) return -1;
        else if (id == irw1) return -1;
        else if (id == irx1) return -1;
        else if (id == iry1) return -1;
        else {
            // An invalid rate id was passed
            AMREX_ASSERT(false);
            return -99;
        }

    }
}


namespace Rates
{
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    static constexpr int lookup_screen_args(const int id)
    {
        // Returns the screening factor index for this rate id
        if (id == ir3a || id == irg3a) {
            return 0;
        } else if (id == ir3a_stage2 || id == irg3a_stage2) {
            return 1;
        } else if (id == ircag || id == iroga) {
            return 2;
        } else if (id == ir1212) {
            return 3;
        } else if (id == ir1216) {
            return 4;
        } else if (id == ir1616) {
            return 5;
        } else if (id == iroag || id == irnega) {
            return 6;
        } else if (id == irneag || id == irmgga) {
            return 7;
        } else if (id == irmgag || id == irsiga || id == irmgap || id == iralpa) {
            return 8;
        } else if (id == iralpg || id == irsigp) {
            return 9;
        } else if (id == irsiag || id == irsga || id == irsiap || id == irppa) {
            return 10;
        } else if (id == irppg || id == irsgp) {
            return 11;
        } else if (id == irsag || id == irarga || id == irsap || id == irclpa) {
            return 12;
        } else if (id == irclpg || id == irargp) {
            return 13;
        } else if (id == irarag || id == ircaga || id == irarap || id == irkpa) {
            return 14;
        } else if (id == irkpg || id == ircagp) {
            return 15;
        } else if (id == ircaag || id == irtiga || id == ircaap || id == irscpa) {
            return 16;
        } else if (id == irscpg || id == irtigp) {
            return 17;
        } else if (id == irtiag || id == ircrga || id == irtiap || id == irvpa) {
            return 18;
        } else if (id == irvpg || id == ircrgp) {
            return 19;
        } else if (id == ircrag || id == irfega || id == ircrap || id == irmnpa) {
            return 20;
        } else if (id == irmnpg || id == irfegp) {
            return 21;
        } else if (id == irfeag || id == irniga || id == irfeap || id == ircopa) {
            return 22;
        } else if (id == ircopg || id == irnigp) {
            return 23;
        } else {
            // if we are trying to look up an invalid screening index
            // then that is a code bug
            AMREX_ASSERT(false);
            return -1;
        }
    }


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real aprox13tab(const Real& btemp, const Real& bden,
                    const int group, const int id)
    {    
        using namespace RateTable;

        // hash locate
        constexpr int mp = 4;
        const int iat1 = static_cast<int>((std::log10(btemp) - tab_tlo)/tab_tstp) + 1;
        const int iat = amrex::max(1, amrex::min(iat - 1, tab_imax - mp + 1));

        // setup the lagrange interpolation coefficients for a cubic
        const Real x  = btemp;
        const Real x1 = ttab(iat);
        const Real x2 = ttab(iat+1);
        const Real x3 = ttab(iat+2);
        const Real x4 = ttab(iat+3);
        const Real a  = x - x1;
        const Real b  = x - x2;
        const Real c  = x - x3;
        const Real d  = x - x4;
        const Real e  = x1 - x2;
        const Real f  = x1 - x3;
        const Real g  = x1 - x4;
        const Real h  = x2 - x3;
        const Real p  = x2 - x4;
        const Real q  = x3 - x4;
        const Real alfa =  b*c*d/(e*f*g);
        const Real beta = -a*c*d/(e*h*p);
        const Real gama =  a*b*d/(f*h*q);
        const Real delt = -a*b*c/(g*p*q);

        const int dpow = density_power(id);

        auto dfac = [](const int p, const Real x) -> Real {
            if (p == 0)
                return 1.0_rt;
            else if (p == 1)
                return x;
            else if (p == 2)
                return x * x;
            else return 0.0_rt;
        };

        const Real density_factor = dfac(dpow, bden);

        if (group == 1) {
        
            return  (  alfa * rattab(id,iat  )
                     + beta * rattab(id,iat+1)
                     + gama * rattab(id,iat+2)
                     + delt * rattab(id,iat+3) ) * density_factor;

        } else if (group == 2) {

            return  (  alfa * drattabdt(id,iat  )
                     + beta * drattabdt(id,iat+1)
                     + gama * drattabdt(id,iat+2)
                     + delt * drattabdt(id,iat+3) ) * density_factor;

        } else {

            return 0.0_rt;

        }
    }


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    Real aprox13rat(const tf_t& tf, const Real& bden,
                    const int group, const int id)
    {
        using namespace Rates;

        // this routine generates unscreened
        // nuclear reaction rates for the aprox13 network.

        if (tf.temp < 1.0e6_rt) return 0.0_rt;

#define FR_RATE(fid, rid, rfun) ({                   \
    if (id == fid || id == rid) {                    \
        Real frate, dfratedt, rrate, drratedt;       \
                                                     \
        rfun(tf,bden,frate,dfratedt,rrate,drratedt); \
                                                     \
        if (id == fid) {                             \
            if (group == 1) return frate;            \
            else if (group == 2) return dfratedt;    \
            else return frate;                       \
        } else if (id == rid) {                      \
            if (group == 1) return rrate;            \
            else if (group == 2) return drratedt;    \
            else return rrate;                       \
        } else return frate;                         \
    }                                                \
})

#define F_RATE(fid, rfun) ({                         \
    if (id == fid) {                                 \
        Real frate, dfratedt, rrate, drratedt;       \
                                                     \
        rfun(tf,bden,frate,dfratedt,rrate,drratedt); \
                                                     \
        if (group == 1) return frate;                \
        else return dfratedt;                        \
    }                                                \
})

        // Determine which c12(a,g)o16 rate to use
        if (use_c12ag_deboer17) {
            // deboer + 2017 c12(a,g)o16 rate
            FR_RATE(ircag, iroga, rate_c12ag_deboer17);
        } else {
            // 1.7 times cf88 c12(a,g)o16 rate
            FR_RATE(ircag, iroga, rate_c12ag);
        }

        // triple alpha to c12
        FR_RATE(ir3a, irg3a, rate_triplealf);

        // c12 + c12
        F_RATE(ir1212, rate_c12c12);

        // c12 + o16
        F_RATE(ir1216, rate_c12o16);

        // o16 + o16
        F_RATE(ir1616, rate_o16o16);

        // o16(a,g)ne20
        FR_RATE(iroag, irnega, rate_o16ag);

        // ne20(a,g)mg24
        FR_RATE(irneag, irmgga, rate_ne20ag);

        // mg24(a,g)si28
        FR_RATE(irmgag, irsiga, rate_mg24ag);

        // mg24(a,p)al27
        FR_RATE(irmgap, iralpa, rate_mg24ap);

        // al27(p,g)si28
        FR_RATE(iralpg, irsigp, rate_al27pg);

        // si28(a,g)s32
        FR_RATE(irsiag, irsga, rate_si28ag);

        // si28(a,p)p31
        FR_RATE(irsiap, irppa, rate_si28ap);

        // p31(p,g)s32
        FR_RATE(irppg, irsgp, rate_p31pg);

        // s32(a,g)ar36
        FR_RATE(irsag, irarga, rate_s32ag);

        // s32(a,p)cl35
        FR_RATE(irsap, irclpa, rate_s32ap);

        // cl35(p,g)ar36
        FR_RATE(irclpg, irargp, rate_cl35pg);

        // ar36(a,g)ca40
        FR_RATE(irarag, ircaga, rate_ar36ag);

        // ar36(a,p)k39
        FR_RATE(irarap, irkpa, rate_ar36ap);

        // k39(p,g)ca40
        FR_RATE(irkpg, ircagp, rate_k39pg);

        // ca40(a,g)ti44
        FR_RATE(ircaag, irtiga, rate_ca40ag);

        // ca40(a,p)sc43
        FR_RATE(ircaap, irscpa, rate_ca40ap);

        // sc43(p,g)ti44
        FR_RATE(irscpg, irtigp, rate_sc43pg);

        // ti44(a,g)cr48
        FR_RATE(irtiag, ircrga, rate_ti44ag);

        // ti44(a,p)v47
        FR_RATE(irtiap, irvpa, rate_ti44ap);

        // v47(p,g)cr48
        FR_RATE(irvpg, ircrgp, rate_v47pg);

        // cr48(a,g)fe52
        FR_RATE(ircrag, irfega, rate_cr48ag);

        // cr48(a,p)mn51
        FR_RATE(ircrap, irmnpa, rate_cr48ap);

        // mn51(p,g)fe52
        FR_RATE(irmnpg, irfegp, rate_mn51pg);

        // fe52(a,g)ni56
        FR_RATE(irfeag, irniga, rate_fe52ag);

        // fe52(a,p)co55
        FR_RATE(irfeap, ircopa, rate_fe52ap);

        // co55(p,g)ni56
        FR_RATE(ircopg, irnigp, rate_co55pg);

        // We should have returned by now if the rate index passed
        // was a valid rate in the network
        AMREX_ASSERT(false);
        return 0.0_rt;
    }

    template<bool use_rate_table, bool use_screening>
    struct RateEvaluatorType
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        RateEvaluatorType(const Real rho, const Real T,
                          const Array1D<Real, 1, NumSpec>& y) :
                          dens(rho), temp(T), tf(get_tfactors(temp)),
                          state(construct_plasma_state(temp, dens, y))
                          {}


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real raw_rates(const int group, const int id) const
        {

            // return the rate value given this group and rate id
            if (use_rate_table) {
                return aprox13tab(temp, dens, group, id);
            } else {
                return aprox13rat(tf, dens, group, id);
            }

        }


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void screen_rates(Real& scf, Real& scfdt,
                        const int group, const int id) const
        {

            if (id == ir3a || id == irg3a) {

                // this is a special case for the 3-body triple alpha reaction
                const int scn_stage1 = lookup_screen_args(id);
                const int scn_stage2 = lookup_screen_args(ir3a_stage2);

                Real sc1a,sc1adt,sc2a,sc2adt,sc1add,sc2add;

                screen(state, scn_stage1, sc1a, sc1adt, sc1add);
                screen(state, scn_stage2, sc2a, sc2adt, sc2add);

                scf   = sc1a * sc2a;
                scfdt = sc1adt*sc2a + sc1a*sc2adt;

            } else {

                Real scfdd;

                const int scn = lookup_screen_args(id);

                screen(state, scn, scf, scfdt, scfdd);

            }

        }


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real primary_rate(const int group, const int id) const
        {

            // primary rates are the non-proton-link rates in the network
            // i.e. the rates that are not approximated

            // first we need the screening factors for this rate id
            Real scf = 1.0_rt;
            Real scfdt = 0.0_rt;

            if (use_screening) {
                screen_rates(scf, scfdt, group, id);
            }

            // return either the screened rate or its temperature derivative
            if (group == 1) {
                // this is the raw rate times the screening factor 
                return raw_rates(group, id) * scf;
            } else {
                // this is d/dT (raw rate * screening factor)
                return raw_rates(group, id) * scf + raw_rates(1, id) * scfdt;  
            }

        }


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real proton_link_rate(const int group, const int id) const
        {

#define PROTON_LINK_RATE(id_link, id_pa, id_pg) ({                        \
    if (id == id_link) {                                                  \
        const Real denom = primary_rate(1,id_pa) + primary_rate(1,id_pg);               \
                                                                          \
        if (denom > 1.0e-30_rt) {                                         \
            const Real zz = 1.0e0_rt/denom;                               \
            const Real prate = primary_rate(1,id_pa)*zz;                         \
            if (group == 1) {                                             \
                return prate;                                             \
            } else {                                                      \
                const Real denomdt = primary_rate(2,id_pa) + primary_rate(2,id_pg);     \
                const Real pratedt = (primary_rate(2,id_pa) - prate*denomdt)*zz; \
                return pratedt;                                           \
            }                                                             \
        } else {                                                          \
            return 0.0_rt;                                                \
        }                                                                 \
    }                                                                     \
})

            // mg24(a,p)27al(p,g)28si
            PROTON_LINK_RATE(irr1, iralpa, iralpg);

            // si28(a,p)p31(p,g)s32
            PROTON_LINK_RATE(irs1, irppa, irppg);

            // s32(a,p)cl35(p,g)ar36
            PROTON_LINK_RATE(irt1, irclpa, irclpg);

            // ar36(a,p)k39(p,g)ca40
            PROTON_LINK_RATE(iru1, irkpa, irkpg);

            // ca40(a,p)sc43(p,g)ti44
            PROTON_LINK_RATE(irv1, irscpa, irscpg);

            // ti44(a,p)v47(p,g)cr48
            PROTON_LINK_RATE(irw1, irvpa, irvpg);

            // cr48(a,p)mn51(p,g)fe52
            PROTON_LINK_RATE(irx1, irmnpa, irmnpg);

            // fe52(a,p)co55(p,g)ni56
            PROTON_LINK_RATE(iry1, ircopa, ircopg);

            // We should have returned by now if the rate
            // index passed was a proton link rate in the network
            AMREX_ASSERT(false);
            return 0.0_rt;
        }


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Real rates(const int group, const int id) const
        {

            // check if we're looking for one of the proton link rates
            if (id >= irr1 && id <= iry1) {
                return proton_link_rate(group, id);
            } else {
                // otherwise, evaluate a non-proton link rate
                return primary_rate(group, id);
            }

        }


        private:
            // thermodynamic conditions at which to evaluate rates
            const Real dens;
            const Real temp;

            // temperature factors for evaluating rates
            const tf_t tf;

            // plasma state for screening rates
            const plasma_state_t state;
    };
}

// typedef depending on whether we are using tables
// (we always have screening turned on)
#ifdef RATE_TABLES
typedef Rates::RateEvaluatorType<true, true> RateEvaluator;
#else
typedef Rates::RateEvaluatorType<false, true> RateEvaluator;
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void set_aprox13rat()
{
    using namespace RateTable;
    using namespace Rates;

    Real btemp;
    Real bden = 1.0e0_rt;

    // The tables do not depend on composition
    // because they are unscreened rates.
    // However, we need to initialize y here
    // because a plasma state is constructed
    // when the RateEvaluator is initialized.
    Array1D<Real, 1, NumSpec> y;

    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = 1.0_rt/NumSpec;
    }

    for (int i = 1; i <= tab_imax; ++i) {

       btemp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
       btemp = std::pow(10.0e0_rt, btemp);

       // Create a RateEvaluatorType with tables and screening disabled
       RateEvaluatorType<false, false> rr(bden, btemp, y);

       ttab(i) = btemp;

       for (int j = 1; j <= Rates::NumRates; ++j) {

          rattab(j,i)    = rr.rates(1,j);
          drattabdt(j,i) = rr.rates(2,j);

       }
    }
}


#define GET_ALL_RATES(INDEX_RATE) \
    const Real GETRATE(ir3a) = rr.rates(INDEX_RATE, ir3a); \
    const Real GETRATE(irg3a) = rr.rates(INDEX_RATE, irg3a); \
    const Real GETRATE(ircag) = rr.rates(INDEX_RATE, ircag); \
    const Real GETRATE(iroga) = rr.rates(INDEX_RATE, iroga); \
    const Real GETRATE(ir1212) = rr.rates(INDEX_RATE, ir1212); \
    const Real GETRATE(ir1216) = rr.rates(INDEX_RATE, ir1216); \
    const Real GETRATE(ir1616) = rr.rates(INDEX_RATE, ir1616); \
    const Real GETRATE(iroag) = rr.rates(INDEX_RATE, iroag); \
    const Real GETRATE(irnega) = rr.rates(INDEX_RATE, irnega); \
    const Real GETRATE(irneag) = rr.rates(INDEX_RATE, irneag); \
    const Real GETRATE(irmgga) = rr.rates(INDEX_RATE, irmgga); \
    const Real GETRATE(irmgag) = rr.rates(INDEX_RATE, irmgag); \
    const Real GETRATE(irsiga) = rr.rates(INDEX_RATE, irsiga); \
    const Real GETRATE(irmgap) = rr.rates(INDEX_RATE, irmgap); \
    const Real GETRATE(iralpa) = rr.rates(INDEX_RATE, iralpa); \
    const Real GETRATE(iralpg) = rr.rates(INDEX_RATE, iralpg); \
    const Real GETRATE(irsigp) = rr.rates(INDEX_RATE, irsigp); \
    const Real GETRATE(irsiag) = rr.rates(INDEX_RATE, irsiag); \
    const Real GETRATE(irsga) = rr.rates(INDEX_RATE, irsga); \
    const Real GETRATE(irsiap) = rr.rates(INDEX_RATE, irsiap); \
    const Real GETRATE(irppa) = rr.rates(INDEX_RATE, irppa); \
    const Real GETRATE(irppg) = rr.rates(INDEX_RATE, irppg); \
    const Real GETRATE(irsgp) = rr.rates(INDEX_RATE, irsgp); \
    const Real GETRATE(irsag) = rr.rates(INDEX_RATE, irsag); \
    const Real GETRATE(irarga) = rr.rates(INDEX_RATE, irarga); \
    const Real GETRATE(irsap) = rr.rates(INDEX_RATE, irsap); \
    const Real GETRATE(irclpa) = rr.rates(INDEX_RATE, irclpa); \
    const Real GETRATE(irclpg) = rr.rates(INDEX_RATE, irclpg); \
    const Real GETRATE(irargp) = rr.rates(INDEX_RATE, irargp); \
    const Real GETRATE(irarag) = rr.rates(INDEX_RATE, irarag); \
    const Real GETRATE(ircaga) = rr.rates(INDEX_RATE, ircaga); \
    const Real GETRATE(irarap) = rr.rates(INDEX_RATE, irarap); \
    const Real GETRATE(irkpa) = rr.rates(INDEX_RATE, irkpa); \
    const Real GETRATE(irkpg) = rr.rates(INDEX_RATE, irkpg); \
    const Real GETRATE(ircagp) = rr.rates(INDEX_RATE, ircagp); \
    const Real GETRATE(ircaag) = rr.rates(INDEX_RATE, ircaag); \
    const Real GETRATE(irtiga) = rr.rates(INDEX_RATE, irtiga); \
    const Real GETRATE(ircaap) = rr.rates(INDEX_RATE, ircaap); \
    const Real GETRATE(irscpa) = rr.rates(INDEX_RATE, irscpa); \
    const Real GETRATE(irscpg) = rr.rates(INDEX_RATE, irscpg); \
    const Real GETRATE(irtigp) = rr.rates(INDEX_RATE, irtigp); \
    const Real GETRATE(irtiag) = rr.rates(INDEX_RATE, irtiag); \
    const Real GETRATE(ircrga) = rr.rates(INDEX_RATE, ircrga); \
    const Real GETRATE(irtiap) = rr.rates(INDEX_RATE, irtiap); \
    const Real GETRATE(irvpa) = rr.rates(INDEX_RATE, irvpa); \
    const Real GETRATE(irvpg) = rr.rates(INDEX_RATE, irvpg); \
    const Real GETRATE(ircrgp) = rr.rates(INDEX_RATE, ircrgp); \
    const Real GETRATE(ircrag) = rr.rates(INDEX_RATE, ircrag); \
    const Real GETRATE(irfega) = rr.rates(INDEX_RATE, irfega); \
    const Real GETRATE(ircrap) = rr.rates(INDEX_RATE, ircrap); \
    const Real GETRATE(irmnpa) = rr.rates(INDEX_RATE, irmnpa); \
    const Real GETRATE(irmnpg) = rr.rates(INDEX_RATE, irmnpg); \
    const Real GETRATE(irfegp) = rr.rates(INDEX_RATE, irfegp); \
    const Real GETRATE(irfeag) = rr.rates(INDEX_RATE, irfeag); \
    const Real GETRATE(irniga) = rr.rates(INDEX_RATE, irniga); \
    const Real GETRATE(irfeap) = rr.rates(INDEX_RATE, irfeap); \
    const Real GETRATE(ircopa) = rr.rates(INDEX_RATE, ircopa); \
    const Real GETRATE(ircopg) = rr.rates(INDEX_RATE, ircopg); \
    const Real GETRATE(irnigp) = rr.rates(INDEX_RATE, irnigp); \
    const Real GETRATE(irr1) = rr.rates(INDEX_RATE, irr1); \
    const Real GETRATE(irs1) = rr.rates(INDEX_RATE, irs1); \
    const Real GETRATE(irt1) = rr.rates(INDEX_RATE, irt1); \
    const Real GETRATE(iru1) = rr.rates(INDEX_RATE, iru1); \
    const Real GETRATE(irv1) = rr.rates(INDEX_RATE, irv1); \
    const Real GETRATE(irw1) = rr.rates(INDEX_RATE, irw1); \
    const Real GETRATE(irx1) = rr.rates(INDEX_RATE, irx1); \
    const Real GETRATE(iry1) = rr.rates(INDEX_RATE, iry1);


#define GETRHS(y, dydt, deriva) ({ \
    for (int i = 1; i <= NumSpec; ++i) { \
        dydt(i) = 0.0_rt; \
    } \
 \
    /* he4 heavy ion reactions */ \
    dydt(He4) += esum3([&](int i) -> Real { \
        if (i == 1) return 0.5e0_rt * y(C12) * y(C12) * GETRATE(ir1212); \
        else if (i == 2) return 0.5e0_rt * y(C12) * y(O16) * GETRATE(ir1216); \
        else if (i == 3) return 0.56e0_rt * 0.5e0_rt * y(O16) * y(O16) * GETRATE(ir1616); \
        else return 0.0_rt; \
    }); \
 \
    /* (a,g) and (g,a) reactions */ \
    dydt(He4) += esum12([&](int i) -> Real { \
        if (i == 1)  return -0.5e0_rt * y(He4) * y(He4) * y(He4) * GETRATE(ir3a); \
        else if (i == 2)  return  3.0e0_rt * y(C12) * GETRATE(irg3a); \
        else if (i == 3)  return -y(He4)  * y(C12) * GETRATE(ircag); \
        else if (i == 4)  return  y(O16)  * GETRATE(iroga); \
        else if (i == 5)  return -y(He4)  * y(O16) * GETRATE(iroag); \
        else if (i == 6)  return  y(Ne20) * GETRATE(irnega); \
        else if (i == 7)  return -y(He4)  * y(Ne20) * GETRATE(irneag); \
        else if (i == 8)  return  y(Mg24) * GETRATE(irmgga); \
        else if (i == 9)  return -y(He4)  * y(Mg24)* GETRATE(irmgag); \
        else if (i == 10) return  y(Si28) * GETRATE(irsiga); \
        else if (i == 11) return -y(He4)  * y(Si28)*GETRATE(irsiag); \
        else if (i == 12) return  y(S32)  * GETRATE(irsga); \
        else return 0.0_rt; \
    }); \
 \
    dydt(He4) += esum12([&](int i) -> Real { \
        if (i == 1)  return -y(He4)  * y(S32) * GETRATE(irsag); \
        else if (i == 2)  return  y(Ar36) * GETRATE(irarga); \
        else if (i == 3)  return -y(He4)  * y(Ar36)*GETRATE(irarag); \
        else if (i == 4)  return  y(Ca40) * GETRATE(ircaga); \
        else if (i == 5)  return -y(He4)  * y(Ca40)*GETRATE(ircaag); \
        else if (i == 6)  return  y(Ti44) * GETRATE(irtiga); \
        else if (i == 7)  return -y(He4)  * y(Ti44)*GETRATE(irtiag); \
        else if (i == 8)  return  y(Cr48) * GETRATE(ircrga); \
        else if (i == 9)  return -y(He4)  * y(Cr48)*GETRATE(ircrag); \
        else if (i == 10) return  y(Fe52) * GETRATE(irfega); \
        else if (i == 11) return -y(He4)  * y(Fe52) * GETRATE(irfeag); \
        else if (i == 12) return  y(Ni56) * GETRATE(irniga); \
        else return 0.0_rt; \
    }); \
 \
    /* (a,p)(p,g) and (g,p)(p,a) reactions */ \
 \
    if (!deriva) { \
 \
        dydt(He4) += esum17([&](int i) -> Real { \
            if (i == 1) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*GETRATE(irs1)*GETRATE(ir1616); \
            else if (i == 2) return -y(He4)  * y(Mg24) * GETRATE(irmgap)*(1.0e0_rt-GETRATE(irr1)); \
            else if (i == 3) return  y(Si28) * GETRATE(irsigp) * GETRATE(irr1); \
            else if (i == 4) return -y(He4)  * y(Si28) * GETRATE(irsiap)*(1.0e0_rt-GETRATE(irs1)); \
            else if (i == 5) return  y(S32)  * GETRATE(irsgp) * GETRATE(irs1); \
            else if (i == 6) return -y(He4)  * y(S32) * GETRATE(irsap)*(1.0e0_rt-GETRATE(irt1)); \
            else if (i == 7) return  y(Ar36) * GETRATE(irargp) * GETRATE(irt1); \
            else if (i == 8) return -y(He4)  * y(Ar36) * GETRATE(irarap)*(1.0e0_rt-GETRATE(iru1)); \
            else if (i == 9) return  y(Ca40) * GETRATE(ircagp) * GETRATE(iru1); \
            else if (i == 10) return -y(He4)  * y(Ca40) * GETRATE(ircaap)*(1.0e0_rt-GETRATE(irv1)); \
            else if (i == 11) return  y(Ti44) * GETRATE(irtigp) * GETRATE(irv1); \
            else if (i == 12) return -y(He4)  * y(Ti44) * GETRATE(irtiap)*(1.0e0_rt-GETRATE(irw1)); \
            else if (i == 13) return  y(Cr48) * GETRATE(ircrgp) * GETRATE(irw1); \
            else if (i == 14) return -y(He4)  * y(Cr48) * GETRATE(ircrap)*(1.0e0_rt-GETRATE(irx1)); \
            else if (i == 15) return  y(Fe52) * GETRATE(irfegp) * GETRATE(irx1); \
            else if (i == 16) return -y(He4)  * y(Fe52) * GETRATE(irfeap)*(1.0e0_rt-GETRATE(iry1)); \
            else if (i == 17) return  y(Ni56) * GETRATE(irnigp) * GETRATE(iry1); \
            else return 0.0_rt; \
        }); \
 \
    } else { \
 \
        dydt(He4) += esum10([&](int i) -> Real { \
            if (i == 1) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rate_irs1 * GETRATE(ir1616); \
            else if (i == 2) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * GETRATE(irs1) * rate_ir1616; \
            else if (i == 3) return -y(He4)*y(Mg24) * GETRATE(irmgap)*(1.0e0_rt - rate_irr1); \
            else if (i == 4) return  y(He4)*y(Mg24) * rate_irmgap*GETRATE(irr1); \
            else if (i == 5) return  y(Si28) * rate_irsigp * GETRATE(irr1); \
            else if (i == 6) return  y(Si28) * GETRATE(irsigp) * rate_irr1; \
            else if (i == 7) return -y(He4)*y(Si28) * GETRATE(irsiap)*(1.0e0_rt - rate_irs1); \
            else if (i == 8) return  y(He4)*y(Si28) * rate_irsiap * GETRATE(irs1); \
            else if (i == 9) return  y(S32)  * rate_irsgp * GETRATE(irs1); \
            else if (i == 10) return  y(S32)  * GETRATE(irsgp) * rate_irs1; \
            else return 0.0_rt; \
        }); \
 \
        dydt(He4) += esum12([&](int i) -> Real { \
            if (i == 1) return -y(He4)*y(S32) * GETRATE(irsap)*(1.0e0_rt - rate_irt1); \
            else if (i == 2) return  y(He4)*y(S32) * rate_irsap*GETRATE(irt1); \
            else if (i == 3) return  y(Ar36) * rate_irargp * GETRATE(irt1); \
            else if (i == 4) return  y(Ar36) * GETRATE(irargp) * rate_irt1; \
            else if (i == 5) return -y(He4)*y(Ar36) * GETRATE(irarap)*(1.0e0_rt - rate_iru1); \
            else if (i == 6) return  y(He4)*y(Ar36) * rate_irarap*GETRATE(iru1); \
            else if (i == 7) return  y(Ca40) * rate_ircagp * GETRATE(iru1); \
            else if (i == 8) return  y(Ca40) * GETRATE(ircagp) * rate_iru1; \
            else if (i == 9) return -y(He4)*y(Ca40) * GETRATE(ircaap)*(1.0e0_rt-rate_irv1); \
            else if (i == 10) return  y(He4)*y(Ca40) * rate_ircaap*GETRATE(irv1); \
            else if (i == 11) return  y(Ti44) * rate_irtigp * GETRATE(irv1); \
            else if (i == 12) return  y(Ti44) * GETRATE(irtigp) * rate_irv1; \
            else return 0.0_rt; \
        }); \
 \
        dydt(He4) += esum12([&](int i) -> Real { \
            if (i == 1) return -y(He4)*y(Ti44) * GETRATE(irtiap)*(1.0e0_rt - rate_irw1); \
            else if (i == 2) return  y(He4)*y(Ti44) * rate_irtiap*GETRATE(irw1); \
            else if (i == 3) return  y(Cr48) * rate_ircrgp * GETRATE(irw1); \
            else if (i == 4) return  y(Cr48) * GETRATE(ircrgp) * rate_irw1; \
            else if (i == 5) return -y(He4)*y(Cr48) * GETRATE(ircrap)*(1.0e0_rt - rate_irx1); \
            else if (i == 6) return  y(He4)*y(Cr48) * rate_ircrap*GETRATE(irx1); \
            else if (i == 7) return  y(Fe52) * rate_irfegp * GETRATE(irx1); \
            else if (i == 8) return  y(Fe52) * GETRATE(irfegp) * rate_irx1; \
            else if (i == 9) return -y(He4)*y(Fe52) * GETRATE(irfeap)*(1.0e0_rt - rate_iry1); \
            else if (i == 10) return  y(He4)*y(Fe52) * rate_irfeap*GETRATE(iry1); \
            else if (i == 11) return  y(Ni56) * rate_irnigp * GETRATE(iry1); \
            else if (i == 12) return  y(Ni56) * GETRATE(irnigp) * rate_iry1; \
            else return 0.0_rt; \
        }); \
     \
    } \
 \
    /* c12 reactions */ \
    dydt(C12) += esum6([&](int i) -> Real { \
        if (i == 1) return -y(C12) * y(C12) * GETRATE(ir1212); \
        else if (i == 2) return -y(C12) * y(O16) * GETRATE(ir1216); \
        else if (i == 3) return  y(He4) * y(He4) * y(He4) * GETRATE(ir3a) / 6.0_rt;; \
        else if (i == 4) return -y(C12) * GETRATE(irg3a); \
        else if (i == 5) return -y(C12) * y(He4) * GETRATE(ircag); \
        else if (i == 6) return  y(O16) * GETRATE(iroga); \
        else return 0.0_rt; \
    }); \
 \
 \
    /* o16 reactions */ \
    dydt(O16) += esum6([&](int i) -> Real { \
        if (i == 1) return -y(C12) * y(O16) * GETRATE(ir1216); \
        else if (i == 2) return -y(O16) * y(O16) * GETRATE(ir1616); \
        else if (i == 3) return  y(C12) * y(He4) * GETRATE(ircag); \
        else if (i == 4) return -y(O16) * y(He4) * GETRATE(iroag); \
        else if (i == 5) return -y(O16) * GETRATE(iroga); \
        else if (i == 6) return  y(Ne20) * GETRATE(irnega); \
        else return 0.0_rt; \
    }); \
 \
 \
    /* ne20 reactions */ \
    dydt(Ne20) += esum5([&](int i) -> Real { \
        if (i == 1) return  0.5e0_rt * y(C12) * y(C12) * GETRATE(ir1212); \
        else if (i == 2) return  y(O16) * y(He4) * GETRATE(iroag); \
        else if (i == 3) return -y(Ne20) * y(He4) * GETRATE(irneag); \
        else if (i == 4) return -y(Ne20) * GETRATE(irnega); \
        else if (i == 5) return  y(Mg24) * GETRATE(irmgga); \
        else return 0.0_rt; \
    }); \
 \
 \
    /* mg24 reactions */ \
    dydt(Mg24) += esum5([&](int i) -> Real { \
        if (i == 1) return  0.5e0_rt * y(C12) * y(O16) * GETRATE(ir1216); \
        else if (i == 2) return  y(Ne20) * y(He4) * GETRATE(irneag); \
        else if (i == 3) return -y(Mg24) * y(He4) * GETRATE(irmgag); \
        else if (i == 4) return -y(Mg24) * GETRATE(irmgga); \
        else if (i == 5) return  y(Si28) * GETRATE(irsiga); \
        else return 0.0_rt; \
    }); \
 \
    if (!deriva) { \
 \
        dydt(Mg24) += -y(Mg24) * y(He4) * GETRATE(irmgap)*(1.0e0_rt-GETRATE(irr1)) + \
                        y(Si28) * GETRATE(irr1) * GETRATE(irsigp); \
 \
    } else { \
 \
        dydt(Mg24) += esum4([&](int i) -> Real { \
            if (i == 1) return -y(Mg24)*y(He4) * GETRATE(irmgap)*(1.0e0_rt - rate_irr1); \
            else if (i == 2) return  y(Mg24)*y(He4) * rate_irmgap*GETRATE(irr1); \
            else if (i == 3) return  y(Si28) * rate_irr1 * GETRATE(irsigp); \
            else if (i == 4) return  y(Si28) * GETRATE(irr1) * rate_irsigp; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
 \
    /* si28 reactions */ \
    dydt(Si28) += esum6([&](int i) -> Real { \
        if (i == 1) return  0.5e0_rt * y(C12) * y(O16) * GETRATE(ir1216); \
        else if (i == 2) return  0.56e0_rt * 0.5e0_rt*y(O16) * y(O16) * GETRATE(ir1616); \
        else if (i == 3) return  y(Mg24) * y(He4) * GETRATE(irmgag); \
        else if (i == 4) return -y(Si28) * y(He4) * GETRATE(irsiag); \
        else if (i == 5) return -y(Si28) * GETRATE(irsiga); \
        else if (i == 6) return  y(S32)  * GETRATE(irsga); \
        else return 0.0_rt; \
    }); \
 \
    if (!deriva) { \
 \
        dydt(Si28) += esum5([&](int i) -> Real { \
            if (i == 1) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*GETRATE(irs1)*GETRATE(ir1616); \
            else if (i == 2) return  y(Mg24) * y(He4) * GETRATE(irmgap)*(1.0e0_rt-GETRATE(irr1)); \
            else if (i == 3) return -y(Si28) * GETRATE(irr1) * GETRATE(irsigp); \
            else if (i == 4) return -y(Si28) * y(He4) * GETRATE(irsiap)*(1.0e0_rt-GETRATE(irs1)); \
            else if (i == 5) return  y(S32)  * GETRATE(irs1) * GETRATE(irsgp); \
            else return 0.0_rt; \
        }); \
 \
    } else { \
 \
        dydt(Si28) += esum10([&](int i) -> Real { \
            if (i == 1) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rate_irs1*GETRATE(ir1616); \
            else if (i == 2) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * GETRATE(irs1)*rate_ir1616; \
            else if (i == 3) return  y(Mg24)*y(He4) * GETRATE(irmgap)*(1.0e0_rt - rate_irr1); \
            else if (i == 4) return -y(Mg24)*y(He4) * rate_irmgap*GETRATE(irr1); \
            else if (i == 5) return -y(Si28) * rate_irr1 * GETRATE(irsigp); \
            else if (i == 6) return -y(Si28) * GETRATE(irr1) * rate_irsigp; \
            else if (i == 7) return -y(Si28)*y(He4) * GETRATE(irsiap)*(1.0e0_rt - rate_irs1); \
            else if (i == 8) return  y(Si28)*y(He4) * rate_irsiap*GETRATE(irs1); \
            else if (i == 9) return y(S32) * rate_irs1 * GETRATE(irsgp); \
            else if (i == 10) return y(S32) * GETRATE(irs1) * rate_irsgp; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
    /* s32 reactions */ \
    dydt(S32) += esum5([&](int i) -> Real { \
        if (i == 1) return  0.1e0_rt * 0.5e0_rt*y(O16) * y(O16) * GETRATE(ir1616); \
        else if (i == 2) return  y(Si28) * y(He4) * GETRATE(irsiag); \
        else if (i == 3) return -y(S32) * y(He4) * GETRATE(irsag); \
        else if (i == 4) return -y(S32) * GETRATE(irsga); \
        else if (i == 5) return  y(Ar36) * GETRATE(irarga); \
        else return 0.0_rt; \
    }); \
 \
 \
    if (!deriva) { \
 \
        dydt(S32) += esum5([&](int i) -> Real { \
            if (i == 1) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)* GETRATE(ir1616)*(1.0e0_rt-GETRATE(irs1)); \
            else if (i == 2) return  y(Si28) * y(He4) * GETRATE(irsiap)*(1.0e0_rt-GETRATE(irs1)); \
            else if (i == 3) return -y(S32) * GETRATE(irs1) * GETRATE(irsgp); \
            else if (i == 4) return -y(S32) * y(He4) * GETRATE(irsap)*(1.0e0_rt-GETRATE(irt1)); \
            else if (i == 5) return  y(Ar36) * GETRATE(irt1) * GETRATE(irargp); \
            else return 0.0_rt; \
        }); \
 \
 \
    } else { \
 \
        dydt(S32) += esum10([&](int i) -> Real { \
            if (i == 1) return  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * GETRATE(ir1616)*(1.0e0_rt-rate_irs1); \
            else if (i == 2) return -0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rate_ir1616*GETRATE(irs1); \
            else if (i == 3) return  y(Si28)*y(He4) * GETRATE(irsiap)*(1.0e0_rt-rate_irs1); \
            else if (i == 4) return -y(Si28)*y(He4) * rate_irsiap*GETRATE(irs1); \
            else if (i == 5) return -y(S32) * rate_irs1 * GETRATE(irsgp); \
            else if (i == 6) return -y(S32) * GETRATE(irs1) * rate_irsgp; \
            else if (i == 7) return -y(S32)*y(He4) * GETRATE(irsap)*(1.0e0_rt-rate_irt1); \
            else if (i == 8) return  y(S32)*y(He4) * rate_irsap*GETRATE(irt1); \
            else if (i == 9) return  y(Ar36) * rate_irt1 * GETRATE(irargp); \
            else if (i == 10) return  y(Ar36) * GETRATE(irt1) * rate_irargp; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
 \
    /* ar36 reactions */ \
    dydt(Ar36) += esum4([&](int i) -> Real { \
        if (i == 1) return  y(S32)  * y(He4) * GETRATE(irsag); \
        else if (i == 2) return -y(Ar36) * y(He4) * GETRATE(irarag); \
        else if (i == 3) return -y(Ar36) * GETRATE(irarga); \
        else if (i == 4) return  y(Ca40) * GETRATE(ircaga); \
        else return 0.0_rt; \
    }); \
 \
 \
    if (!deriva) { \
 \
        dydt(Ar36) += esum4([&](int i) -> Real { \
            if (i == 1) return y(S32)  * y(He4) * GETRATE(irsap)*(1.0e0_rt-GETRATE(irt1)); \
            else if (i == 2) return -y(Ar36) * GETRATE(irt1) * GETRATE(irargp); \
            else if (i == 3) return -y(Ar36) * y(He4) * GETRATE(irarap)*(1.0e0_rt-GETRATE(iru1)); \
            else if (i == 4) return  y(Ca40) * GETRATE(ircagp) * GETRATE(iru1); \
            else return 0.0_rt; \
        }); \
 \
    } else { \
 \
        dydt(Ar36) += esum8([&](int i) -> Real { \
            if (i == 1) return  y(S32)*y(He4) * GETRATE(irsap)*(1.0e0_rt - rate_irt1); \
            else if (i == 2) return -y(S32)*y(He4) * rate_irsap*GETRATE(irt1); \
            else if (i == 3) return -y(Ar36) * rate_irt1 * GETRATE(irargp); \
            else if (i == 4) return -y(Ar36) * GETRATE(irt1) * rate_irargp; \
            else if (i == 5) return -y(Ar36)*y(He4) * GETRATE(irarap)*(1.0e0_rt-rate_iru1); \
            else if (i == 6) return  y(Ar36)*y(He4) * rate_irarap*GETRATE(iru1); \
            else if (i == 7) return  y(Ca40) * rate_ircagp * GETRATE(iru1); \
            else if (i == 8) return  y(Ca40) * GETRATE(ircagp) * rate_iru1; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
 \
    /* ca40 reactions */ \
    dydt(Ca40) += esum4([&](int i) -> Real { \
        if (i == 1) return  y(Ar36) * y(He4) * GETRATE(irarag); \
        else if (i == 2) return -y(Ca40) * y(He4) * GETRATE(ircaag); \
        else if (i == 3) return -y(Ca40) * GETRATE(ircaga); \
        else if (i == 4) return  y(Ti44) * GETRATE(irtiga); \
        else return 0.0_rt; \
    }); \
 \
 \
    if (!deriva) { \
 \
        dydt(Ca40) += esum4([&](int i) -> Real { \
            if (i == 1) return  y(Ar36) * y(He4) * GETRATE(irarap)*(1.0e0_rt-GETRATE(iru1)); \
            else if (i == 2) return -y(Ca40) * GETRATE(ircagp) * GETRATE(iru1); \
            else if (i == 3) return -y(Ca40) * y(He4) * GETRATE(ircaap)*(1.0e0_rt-GETRATE(irv1)); \
            else if (i == 4) return  y(Ti44) * GETRATE(irtigp) * GETRATE(irv1); \
            else return 0.0_rt; \
        }); \
 \
    } else { \
 \
        dydt(Ca40) += esum8([&](int i) -> Real { \
            if (i == 1) return  y(Ar36)*y(He4) * GETRATE(irarap)*(1.0e0_rt-rate_iru1); \
            else if (i == 2) return -y(Ar36)*y(He4) * rate_irarap*GETRATE(iru1); \
            else if (i == 3) return -y(Ca40) * rate_ircagp * GETRATE(iru1); \
            else if (i == 4) return -y(Ca40) * GETRATE(ircagp) * rate_iru1; \
            else if (i == 5) return -y(Ca40)*y(He4) * GETRATE(ircaap)*(1.0e0_rt-rate_irv1); \
            else if (i == 6) return  y(Ca40)*y(He4) * rate_ircaap*GETRATE(irv1); \
            else if (i == 7) return  y(Ti44) * rate_irtigp * GETRATE(irv1); \
            else if (i == 8) return  y(Ti44) * GETRATE(irtigp) * rate_irv1; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
 \
    /* ti44 reactions */ \
    dydt(Ti44) += esum4([&](int i) -> Real { \
        if (i == 1) return  y(Ca40) * y(He4) * GETRATE(ircaag); \
        else if (i == 2) return -y(Ti44) * y(He4) * GETRATE(irtiag); \
        else if (i == 3) return -y(Ti44) * GETRATE(irtiga); \
        else if (i == 4) return  y(Cr48) * GETRATE(ircrga); \
        else return 0.0_rt; \
    }); \
 \
 \
    if (!deriva) { \
 \
        dydt(Ti44) += esum4([&](int i) -> Real { \
            if (i == 1) return  y(Ca40) * y(He4) * GETRATE(ircaap)*(1.0e0_rt-GETRATE(irv1)); \
            else if (i == 2) return -y(Ti44) * GETRATE(irv1) * GETRATE(irtigp); \
            else if (i == 3) return -y(Ti44) * y(He4) * GETRATE(irtiap)*(1.0e0_rt-GETRATE(irw1)); \
            else if (i == 4) return  y(Cr48) * GETRATE(irw1) * GETRATE(ircrgp); \
            else return 0.0_rt; \
        }); \
 \
    } else { \
 \
        dydt(Ti44) += esum8([&](int i) -> Real { \
            if (i == 1) return  y(Ca40)*y(He4) * GETRATE(ircaap)*(1.0e0_rt-rate_irv1); \
            else if (i == 2) return -y(Ca40)*y(He4) * rate_ircaap*GETRATE(irv1); \
            else if (i == 3) return -y(Ti44) * rate_irv1 * GETRATE(irtigp); \
            else if (i == 4) return -y(Ti44) * GETRATE(irv1) * rate_irtigp; \
            else if (i == 5) return -y(Ti44)*y(He4) * GETRATE(irtiap)*(1.0e0_rt-rate_irw1); \
            else if (i == 6) return  y(Ti44)*y(He4) * rate_irtiap*GETRATE(irw1); \
            else if (i == 7) return  y(Cr48) * rate_irw1 * GETRATE(ircrgp); \
            else if (i == 8) return  y(Cr48) * GETRATE(irw1) * rate_ircrgp; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
 \
    /* cr48 reactions */ \
    dydt(Cr48) += esum4([&](int i) -> Real { \
        if (i == 1) return  y(Ti44) * y(He4) * GETRATE(irtiag); \
        else if (i == 2) return -y(Cr48) * y(He4) * GETRATE(ircrag); \
        else if (i == 3) return -y(Cr48) * GETRATE(ircrga); \
        else if (i == 4) return  y(Fe52) * GETRATE(irfega); \
        else return 0.0_rt; \
    }); \
 \
    if (!deriva) { \
 \
        dydt(Cr48) += esum4([&](int i) -> Real { \
            if (i == 1) return  y(Ti44) * y(He4) * GETRATE(irtiap)*(1.0e0_rt-GETRATE(irw1)); \
            else if (i == 2) return -y(Cr48) * GETRATE(irw1) * GETRATE(ircrgp); \
            else if (i == 3) return -y(Cr48) * y(He4) * GETRATE(ircrap)*(1.0e0_rt-GETRATE(irx1)); \
            else if (i == 4) return  y(Fe52) * GETRATE(irx1) * GETRATE(irfegp); \
            else return 0.0_rt; \
        }); \
 \
    } else { \
 \
        dydt(Cr48) += esum8([&](int i) -> Real { \
            if (i == 1) return  y(Ti44)*y(He4) * GETRATE(irtiap)*(1.0e0_rt-rate_irw1); \
            else if (i == 2) return -y(Ti44)*y(He4) * rate_irtiap*GETRATE(irw1); \
            else if (i == 3) return -y(Cr48) * rate_irw1 * GETRATE(ircrgp); \
            else if (i == 4) return -y(Cr48) * GETRATE(irw1) * rate_ircrgp; \
            else if (i == 5) return -y(Cr48)*y(He4) * GETRATE(ircrap)*(1.0e0_rt-rate_irx1); \
            else if (i == 6) return  y(Cr48)*y(He4) * rate_ircrap*GETRATE(irx1); \
            else if (i == 7) return  y(Fe52) * rate_irx1 * GETRATE(irfegp); \
            else if (i == 8) return  y(Fe52) * GETRATE(irx1) * rate_irfegp; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
    /* fe52 reactions */ \
    dydt(Fe52) += esum4([&](int i) -> Real { \
        if (i == 1) return  y(Cr48) * y(He4) * GETRATE(ircrag); \
        else if (i == 2) return -y(Fe52) * y(He4) * GETRATE(irfeag); \
        else if (i == 3) return -y(Fe52) * GETRATE(irfega); \
        else if (i == 4) return  y(Ni56) * GETRATE(irniga); \
        else return 0.0_rt; \
    }); \
 \
    if (!deriva) { \
 \
        dydt(Fe52) += esum4([&](int i) -> Real { \
            if (i == 1) return  y(Cr48) * y(He4) * GETRATE(ircrap)*(1.0e0_rt-GETRATE(irx1)); \
            else if (i == 2) return -y(Fe52) * GETRATE(irx1) * GETRATE(irfegp); \
            else if (i == 3) return -y(Fe52) * y(He4) * GETRATE(irfeap)*(1.0e0_rt-GETRATE(iry1)); \
            else if (i == 4) return  y(Ni56) * GETRATE(iry1) * GETRATE(irnigp); \
            else return 0.0_rt; \
        }); \
 \
    } else { \
 \
        dydt(Fe52) += esum8([&](int i) -> Real { \
            if (i == 1) return  y(Cr48)*y(He4) * GETRATE(ircrap)*(1.0e0_rt-rate_irx1); \
            else if (i == 2) return -y(Cr48)*y(He4) * rate_ircrap*GETRATE(irx1); \
            else if (i == 3) return -y(Fe52) * rate_irx1 * GETRATE(irfegp); \
            else if (i == 4) return -y(Fe52) * GETRATE(irx1) * rate_irfegp; \
            else if (i == 5) return -y(Fe52)*y(He4) * GETRATE(irfeap)*(1.0e0_rt-rate_iry1); \
            else if (i == 6) return  y(Fe52)*y(He4) * rate_irfeap*GETRATE(iry1); \
            else if (i == 7) return  y(Ni56) * rate_iry1 * GETRATE(irnigp); \
            else if (i == 8) return  y(Ni56) * GETRATE(iry1) * rate_irnigp; \
            else return 0.0_rt; \
        }); \
 \
    } \
 \
    /* ni56 reactions */ \
    dydt(Ni56) += y(Fe52) * y(He4) * GETRATE(irfeag) + \
                 -y(Ni56) * GETRATE(irniga); \
 \
    if (!deriva) { \
 \
        dydt(Ni56) += y(Fe52) * y(He4) * GETRATE(irfeap)*(1.0e0_rt-GETRATE(iry1)) + \
                     -y(Ni56) * GETRATE(iry1) * GETRATE(irnigp); \
 \
    } else { \
 \
        dydt(Ni56) += esum4([&](int i) -> Real { \
            if (i == 1) return  y(Fe52)*y(He4) * GETRATE(irfeap)*(1.0e0_rt-rate_iry1); \
            else if (i == 2) return -y(Fe52)*y(He4) * rate_irfeap*GETRATE(iry1); \
            else if (i == 3) return -y(Ni56) * rate_iry1 * GETRATE(irnigp); \
            else if (i == 4) return -y(Ni56) * GETRATE(iry1) * rate_irnigp; \
            else return 0.0_rt; \
        }); \
 \
    } \
})


// Evaluates the species right hand side of the aprox13 ODEs
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void Fy_aprox13(Array1D<Real, 1, NumSpec> const& y,
                RateEvaluator const& rr,
                Array1D<Real, 1, neqs>& dydt)
{
    using namespace Species;
    using namespace Rates;

#undef GETRATE
#define GETRATE(b) rate_ ## b
    GET_ALL_RATES(1)

    constexpr bool deriva = false;

    GETRHS(y, dydt, deriva);
#undef GETRATE

}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dFydy_dFydT_aprox13(Array1D<Real, 1, NumSpec> const& y,
                         RateEvaluator const& rr,
                         MatrixType& jac)
{

    // this routine sets up the aprox13 jacobian for the isotopes

    using namespace Species;
    using namespace Rates;

    // declare both rates and rate temperature derivatives as const Reals
#undef GETRATE
#define GETRATE(b) rate_ ## b
GET_ALL_RATES(1)
#undef GETRATE
#define GETRATE(b) drate_dT_ ## b
GET_ALL_RATES(2)
#undef GETRATE

    // d(he4)/d(he4)
    jac(He4,He4) = esum20([&](int i) -> Real {
        if (i == 1) return -1.5e0_rt * y(He4) * y(He4) * rate_ir3a;
        else if (i == 2) return -y(C12)  * rate_ircag;
        else if (i == 3) return -y(O16)  * rate_iroag;
        else if (i == 4) return -y(Ne20) * rate_irneag;
        else if (i == 5) return -y(Mg24) * rate_irmgag;
        else if (i == 6) return -y(Si28) * rate_irsiag;
        else if (i == 7) return -y(S32)  * rate_irsag;
        else if (i == 8) return -y(Ar36) * rate_irarag;
        else if (i == 9) return -y(Ca40) * rate_ircaag;
        else if (i == 10) return -y(Ti44) * rate_irtiag;
        else if (i == 11) return -y(Cr48) * rate_ircrag;
        else if (i == 12) return -y(Fe52) * rate_irfeag;
        else if (i == 13) return -y(Mg24) * rate_irmgap * (1.0e0_rt-rate_irr1);
        else if (i == 14) return -y(Si28) * rate_irsiap * (1.0e0_rt-rate_irs1);
        else if (i == 15) return -y(S32)  * rate_irsap  * (1.0e0_rt-rate_irt1);
        else if (i == 16) return -y(Ar36) * rate_irarap * (1.0e0_rt-rate_iru1);
        else if (i == 17) return -y(Ca40) * rate_ircaap * (1.0e0_rt-rate_irv1);
        else if (i == 18) return -y(Ti44) * rate_irtiap * (1.0e0_rt-rate_irw1);
        else if (i == 19) return -y(Cr48) * rate_ircrap * (1.0e0_rt-rate_irx1);
        else if (i == 20) return -y(Fe52) * rate_irfeap * (1.0e0_rt-rate_iry1);
        else return 0.0_rt;
    });

    // (he4)/d(o16)
    jac(He4,O16) = esum5([&](int i) -> Real {
        if (i == 1) return  0.5e0_rt * y(C12) * rate_ir1216;
        else if (i == 2) return  1.12e0_rt * 0.5e0_rt*y(O16) * rate_ir1616;
        else if (i == 3) return  0.68e0_rt * rate_irs1 * 0.5e0_rt*y(O16) * rate_ir1616;
        else if (i == 4) return  rate_iroga;
        else if (i == 5) return -y(He4) * rate_iroag;
        else return 0.0_rt;
    });

    // d(he4)/d(c12)
    jac(He4,C12) = esum4([&](int i) -> Real {
        if (i == 1) return  y(C12) * rate_ir1212;
        else if (i == 2) return  0.5e0_rt * y(O16) * rate_ir1216;
        else if (i == 3) return  3.0e0_rt * rate_irg3a;
        else if (i == 4) return -y(He4) * rate_ircag;
        else return 0.0_rt;
    });

    // d(he4)/d(si28)
    jac(He4,Si28) = esum4([&](int i) -> Real {
        if (i == 1) return  rate_irsiga;
        else if (i == 2) return -y(He4) * rate_irsiag;
        else if (i == 3) return -y(He4) * rate_irsiap * (1.0e0_rt-rate_irs1);
        else if (i == 4) return  rate_irr1 * rate_irsigp;
        else return 0.0_rt;
    });

    // d(he4)/d(s32)
    jac(He4,S32) = esum4([&](int i) -> Real {
        if (i == 1) return  rate_irsga;
        else if (i == 2) return -y(He4) * rate_irsag;
        else if (i == 3) return -y(He4) * rate_irsap * (1.0e0_rt-rate_irt1);
        else if (i == 4) return  rate_irs1 * rate_irsgp;
        else return 0.0_rt;
    });

    // d(he4)/d(ar36)
    jac(He4,Ar36) = esum4([&](int i) -> Real {
        if (i == 1) return  rate_irarga;
        else if (i == 2) return -y(He4) * rate_irarag;
        else if (i == 3) return -y(He4) * rate_irarap * (1.0e0_rt-rate_iru1);
        else if (i == 4) return  rate_irt1 * rate_irargp;
        else return 0.0_rt;
    });

    // d(he4)/d(ca40)
    jac(He4,Ca40) = esum4([&](int i) -> Real {
        if (i == 1) return  rate_ircaga;
        else if (i == 2) return -y(He4) * rate_ircaag;
        else if (i == 3) return -y(He4) * rate_ircaap * (1.0e0_rt-rate_irv1);
        else if (i == 4) return  rate_iru1 * rate_ircagp;
        else return 0.0_rt;
    });

    // d(he4)/d(ti44)
    jac(He4,Ti44) = esum4([&](int i) -> Real {
        if (i == 1) return  rate_irtiga;
        else if (i == 2) return -y(He4) * rate_irtiag;
        else if (i == 3) return -y(He4) * rate_irtiap * (1.0e0_rt-rate_irw1);
        else if (i == 4) return  rate_irv1 * rate_irtigp;
        else return 0.0_rt;
    });

    // d(he4)/d(cr48)
    jac(He4,Cr48) = esum4([&](int i) -> Real {
        if (i == 1) return  rate_ircrga;
        else if (i == 2) return -y(He4) * rate_ircrag;
        else if (i == 3) return -y(He4) * rate_ircrap * (1.0e0_rt-rate_irx1);
        else if (i == 4) return  rate_irw1 * rate_ircrgp;
        else return 0.0_rt;
    });

    // d(he4)/d(fe52)
    jac(He4,Fe52) = esum4([&](int i) -> Real {
        if (i == 1) return  rate_irfega;
        else if (i == 2) return -y(He4) * rate_irfeag;
        else if (i == 3) return -y(He4) * rate_irfeap * (1.0e0_rt-rate_iry1);
        else if (i == 4) return  rate_irx1 * rate_irfegp;
        else return 0.0_rt;
    });

    // d(c12)/d(c12)
    jac(C12,C12) = esum4([&](int i) -> Real {
        if (i == 1) return -2.0e0_rt * y(C12) * rate_ir1212;
        else if (i == 2) return -y(O16) * rate_ir1216;
        else if (i == 3) return -rate_irg3a;
        else if (i == 4) return -y(He4) * rate_ircag;
        else return 0.0_rt;
    });

    // d(o16)/d(o16)
    jac(O16,O16) = esum4([&](int i) -> Real {
        if (i == 1) return -y(C12) * rate_ir1216;
        else if (i == 2) return -2.0e0_rt * y(O16) * rate_ir1616;
        else if (i == 3) return -y(He4) * rate_iroag;
        else if (i == 4) return -rate_iroga;
        else return 0.0_rt;
    });

    // d(si28)/d(he4)
    jac(Si28,He4) = esum4([&](int i) -> Real {
        if (i == 1) return  y(Mg24) * rate_irmgag;
        else if (i == 2) return -y(Si28) * rate_irsiag;
        else if (i == 3) return  y(Mg24) * rate_irmgap * (1.0e0_rt-rate_irr1);
        else if (i == 4) return -y(Si28) * rate_irsiap * (1.0e0_rt-rate_irs1);
        else return 0.0_rt;
    });

    // d(si28)/d(si28)
    jac(Si28,Si28) = esum4([&](int i) -> Real {
        if (i == 1) return  -y(He4) * rate_irsiag;
        else if (i == 2) return -rate_irsiga;
        else if (i == 3) return -rate_irr1 * rate_irsigp;
        else if (i == 4) return -y(He4) * rate_irsiap * (1.0e0_rt-rate_irs1);
        else return 0.0_rt;
    });

    // d(s32)/d(he4)
    jac(S32,He4) = esum4([&](int i) -> Real {
        if (i == 1) return  y(Si28) * rate_irsiag;
        else if (i == 2) return -y(S32) * rate_irsag;
        else if (i == 3) return  y(Si28) * rate_irsiap * (1.0e0_rt-rate_irs1);
        else if (i == 4) return -y(S32) * rate_irsap * (1.0e0_rt-rate_irt1);
        else return 0.0_rt;
    });

    // d(s32)/d(s32)
    jac(S32,S32) = esum4([&](int i) -> Real {
        if (i == 1) return -y(He4) * rate_irsag;
        else if (i == 2) return -rate_irsga;
        else if (i == 3) return -rate_irs1 * rate_irsgp;
        else if (i == 4) return -y(He4) * rate_irsap * (1.0e0_rt-rate_irt1);
        else return 0.0_rt;
    });

    // d(ar36)/d(he4)
    jac(Ar36,He4) = esum4([&](int i) -> Real {
        if (i == 1) return  y(S32)  * rate_irsag;
        else if (i == 2) return -y(Ar36) * rate_irarag;
        else if (i == 3) return  y(S32)  * rate_irsap * (1.0e0_rt-rate_irt1);
        else if (i == 4) return -y(Ar36) * rate_irarap * (1.0e0_rt-rate_iru1);
        else return 0.0_rt;
    });

    // d(ar36)/d(ar36)
    jac(Ar36,Ar36) = esum4([&](int i) -> Real {
        if (i == 1) return -y(He4) * rate_irarag;
        else if (i == 2) return -rate_irarga;
        else if (i == 3) return -rate_irt1 * rate_irargp;
        else if (i == 4) return -y(He4) * rate_irarap * (1.0e0_rt-rate_iru1);
        else return 0.0_rt;
    });

    // d(ca40)/d(he4)
    jac(Ca40,He4) = esum4([&](int i) -> Real {
        if (i == 1) return  y(Ar36) * rate_irarag;
        else if (i == 2) return -y(Ca40) * rate_ircaag;
        else if (i == 3) return  y(Ar36) * rate_irarap*(1.0e0_rt-rate_iru1);
        else if (i == 4) return -y(Ca40) * rate_ircaap*(1.0e0_rt-rate_irv1);
        else return 0.0_rt;
    });

    // d(ca40)/d(ca40)
    jac(Ca40,Ca40) = esum4([&](int i) -> Real {
        if (i == 1) return -y(He4) * rate_ircaag;
        else if (i == 2) return -rate_ircaga;
        else if (i == 3) return -rate_ircagp * rate_iru1;
        else if (i == 4) return -y(He4) * rate_ircaap*(1.0e0_rt-rate_irv1);
        else return 0.0_rt;
    });

    // d(ti44)/d(he4)
    jac(Ti44,He4) = esum4([&](int i) -> Real {
        if (i == 1) return  y(Ca40) * rate_ircaag;
        else if (i == 2) return -y(Ti44) * rate_irtiag;
        else if (i == 3) return  y(Ca40) * rate_ircaap*(1.0e0_rt-rate_irv1);
        else if (i == 4) return -y(Ti44) * rate_irtiap*(1.0e0_rt-rate_irw1);
        else return 0.0_rt;
    });

    // d(ti44)/d(ti44)
    jac(Ti44,Ti44) = esum4([&](int i) -> Real {
        if (i == 1) return -y(He4) * rate_irtiag;
        else if (i == 2) return -rate_irtiga;
        else if (i == 3) return -rate_irv1 * rate_irtigp;
        else if (i == 4) return -y(He4) * rate_irtiap*(1.0e0_rt-rate_irw1);
        else return 0.0_rt;
    });

    // d(cr48)/d(he4)
    jac(Cr48,He4) = esum4([&](int i) -> Real {
        if (i == 1) return  y(Ti44) * rate_irtiag;
        else if (i == 2) return -y(Cr48) * rate_ircrag;
        else if (i == 3) return  y(Ti44) * rate_irtiap*(1.0e0_rt-rate_irw1);
        else if (i == 4) return -y(Cr48) * rate_ircrap*(1.0e0_rt-rate_irx1);
        else return 0.0_rt;
    });

    // d(cr48)/d(cr48)
    jac(Cr48,Cr48) = esum4([&](int i) -> Real {
        if (i == 1) return -y(He4) * rate_ircrag;
        else if (i == 2) return -rate_ircrga;
        else if (i == 3) return -rate_irw1 * rate_ircrgp;
        else if (i == 4) return -y(He4) * rate_ircrap*(1.0e0_rt-rate_irx1);
        else return 0.0_rt;
    });

    // d(fe52)/d(he4)
    jac(Fe52,He4) = esum4([&](int i) -> Real {
        if (i == 1) return  y(Cr48) * rate_ircrag;
        else if (i == 2) return -y(Fe52) * rate_irfeag;
        else if (i == 3) return  y(Cr48) * rate_ircrap * (1.0e0_rt-rate_irx1);
        else if (i == 4) return -y(Fe52) * rate_irfeap * (1.0e0_rt-rate_iry1);
        else return 0.0_rt;
    });

    // d(fe52)/d(fe52)
    jac(Fe52,Fe52) = esum4([&](int i) -> Real {
        if (i == 1) return -y(He4) * rate_irfeag;
        else if (i == 2) return -rate_irfega;
        else if (i == 3) return -rate_irx1 * rate_irfegp;
        else if (i == 4) return -y(He4) * rate_irfeap * (1.0e0_rt-rate_iry1);
        else return 0.0_rt;
    });

    // d(he4)/d(mg24)
    jac(He4,Mg24) = esum3([&](int i) -> Real {
        if (i == 1) return  rate_irmgga;
        else if (i == 2) return -y(He4) * rate_irmgag;
        else if (i == 3) return -y(He4) * rate_irmgap * (1.0e0_rt-rate_irr1);
        else return 0.0_rt;
    });

    // d(mg24)/d(he4)
    jac(Mg24,He4) = esum3([&](int i) -> Real {
        if (i == 1) return  y(Ne20) * rate_irneag;
        else if (i == 2) return -y(Mg24) * rate_irmgag;
        else if (i == 3) return -y(Mg24) * rate_irmgap * (1.0e0_rt-rate_irr1);
        else return 0.0_rt;
    });

    // d(mg24)/d(mg24)
    jac(Mg24,Mg24) = esum3([&](int i) -> Real {
        if (i == 1) return -y(He4) * rate_irmgag;
        else if (i == 2) return -rate_irmgga;
        else if (i == 3) return -y(He4) * rate_irmgap * (1.0e0_rt-rate_irr1);
        else return 0.0_rt;
    });

    // d(si28)/d(o16)
    jac(Si28,O16) = esum3([&](int i) -> Real {
        if (i == 1) return 0.5e0_rt * y(C12) * rate_ir1216;
        else if (i == 2) return 1.12e0_rt * 0.5e0_rt*y(O16) * rate_ir1616;
        else if (i == 3) return 0.68e0_rt * 0.5e0_rt*y(O16) * rate_irs1 * rate_ir1616;
        else return 0.0_rt;
    });

    // d(he4)/d(ne20)
    jac(He4,Ne20)  =  rate_irnega + 
                      -y(He4) * rate_irneag;

    // d(he4)/d(ni56)
    jac(He4,Ni56) = rate_irniga + 
                    rate_iry1 * rate_irnigp;

    // d(c12)/d(he4)
    jac(C12,He4) = 0.5e0_rt * y(He4) * y(He4) * rate_ir3a +
                   -y(C12) * rate_ircag;

    // d(c12)/d(o16)
    jac(C12,O16) = -y(C12) * rate_ir1216 +
                   rate_iroga;

    // d(o16)/d(he4)
    jac(O16,He4) =  y(C12)*rate_ircag +
                    -y(O16)*rate_iroag;

    // d(o16)/d(c12)
    jac(O16,C12) = -y(O16)*rate_ir1216 +
                    y(He4)*rate_ircag;

    // d(ne20)/d(he4)
    jac(Ne20,He4) =  y(O16) * rate_iroag +
                    -y(Ne20) * rate_irneag;

    // d(ne20)/d(ne20)
    jac(Ne20,Ne20) = -y(He4) * rate_irneag +
                    -rate_irnega;

    // d(mg24)/d(si28)
    jac(Mg24,Si28) = rate_irsiga +
                    rate_irr1 * rate_irsigp;

    // d(si28)/d(mg24)
    jac(Si28,Mg24) =  y(He4) * rate_irmgag +
                    y(He4) * rate_irmgap * (1.0e0_rt-rate_irr1);

    // d(si28)/d(s32)
    jac(Si28,S32) = rate_irsga +
                    rate_irs1 * rate_irsgp;

    // d(s32)/d(o16)
    jac(S32,O16) = 0.68e0_rt*0.5e0_rt*y(O16)*rate_ir1616*(1.0e0_rt-rate_irs1) +
                    0.2e0_rt * 0.5e0_rt*y(O16) * rate_ir1616;

    // d(s32)/d(si28)
    jac(S32,Si28) =y(He4) * rate_irsiag +
                    y(He4) * rate_irsiap * (1.0e0_rt-rate_irs1);

    // d(s32)/d(ar36)
    jac(S32,Ar36) = rate_irarga +
                    rate_irt1 * rate_irargp;

    // d(ar36)/d(s32)
    jac(Ar36,S32) = y(He4) * rate_irsag +
                    y(He4) * rate_irsap * (1.0e0_rt-rate_irt1);

    // d(ar36)/d(ca40)
    jac(Ar36,Ca40) = rate_ircaga +
                    rate_ircagp * rate_iru1;

    // d(ca40)/d(ar36)
    jac(Ca40,Ar36) =  y(He4) * rate_irarag +
                    y(He4) * rate_irarap*(1.0e0_rt-rate_iru1);

    // d(ca40)/d(ti44)
    jac(Ca40,Ti44) = rate_irtiga +
                    rate_irtigp * rate_irv1;

    // d(ti44)/d(ca40)
    jac(Ti44,Ca40) =  y(He4) * rate_ircaag +
                    y(He4) * rate_ircaap*(1.0e0_rt-rate_irv1);

    // d(ti44)/d(cr48)
    jac(Ti44,Cr48) = rate_ircrga +
                    rate_irw1 * rate_ircrgp;

    // d(cr48)/d(ti44)
    jac(Cr48,Ti44) =  y(He4) * rate_irtiag +
                    y(He4) * rate_irtiap*(1.0e0_rt-rate_irw1);

    // d(cr48)/d(fe52)
    jac(Cr48,Fe52) = rate_irfega +
                    rate_irx1 * rate_irfegp;

    // d(fe52)/d(cr48)
    jac(Fe52,Cr48) = y(He4) * rate_ircrag +
                    y(He4) * rate_ircrap * (1.0e0_rt-rate_irx1);

    // d(fe52)/d(ni56)
    jac(Fe52,Ni56) = rate_irniga +
                    rate_iry1 * rate_irnigp;

    // d(ni56)/d(he4)
    jac(Ni56,He4) =  y(Fe52) * rate_irfeag +
                    y(Fe52) * rate_irfeap * (1.0e0_rt-rate_iry1);

    // d(ni56)/d(fe52)
    jac(Ni56,Fe52) = y(He4) * rate_irfeag +
                    y(He4) * rate_irfeap * (1.0e0_rt-rate_iry1);

    // d(ni56)/d(ni56)
    jac(Ni56,Ni56) = -rate_irniga +
                    -rate_iry1 * rate_irnigp;

    // d(o16)/d(ne20)
    jac(O16,Ne20) = rate_irnega;

    // d(ne20)/d(c12)
    jac(Ne20,C12) = y(C12) * rate_ir1212;

    // d(ne20)/d(o16)
    jac(Ne20,O16) = y(He4) * rate_iroag;

    // d(ne20)/d(mg24)
    jac(Ne20,Mg24) = rate_irmgga;

    // d(mg24)/d(c12)
    jac(Mg24,C12) = 0.5e0_rt * y(O16) * rate_ir1216;

    // d(mg24)/d(o16)
    jac(Mg24,O16) = 0.5e0_rt * y(C12) * rate_ir1216;

    // d(mg24)/d(ne20)
    jac(Mg24,Ne20) = y(He4) * rate_irneag;

    // d(si28)/d(c12)
    jac(Si28,C12) = 0.5e0_rt * y(O16) * rate_ir1216;

    // now, fill the J(species, temperature) terms in the jacobian
    auto jac_species_temp = [&](int i) -> Real& { 
        return jac(i, net_itemp);
    };

#undef GETRATE
#define GETRATE(b) drate_dT_ ## b
    constexpr bool deriva = true;
    GETRHS(y, jac_species_temp, deriva);
#undef GETRATE

}


template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ener_gener_rate(T const& dydt)
{

    using namespace aprox13;

    // Computes the instantaneous energy generation rate

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives
    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2
    return Xdot * C::Legacy::enuc_conv2;

}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)
{

    /*
     This routine sets up the system of ODE's for the aprox13
     nuclear reactions.  This is an alpha chain + heavy ion network
     with (a,p)(p,g) links.

     Isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
               ar36, ca40, ti44, cr48, fe52, ni56
    */

    // Initialize ydot to 0

    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }

    // Get the data from the state

    const Real rho  = state.rho;
    const Real temp = state.T;
    const Real abar = state.abar;
    const Real zbar = state.zbar;

    Array1D<Real, 1, NumSpec> y;
    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    const RateEvaluator rr(rho, temp, y);

    // Get species RHS

    Fy_aprox13(y, rr, ydot);

    // Instantaneous energy generation rate -- this needs molar fractions

    const Real enuc = ener_gener_rate(ydot);

    // Get the neutrino losses

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Set the energy RHS (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

    // Set the temperature RHS

    temperature_rhs(state, ydot);

}


// Analytical Jacobian
template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void actual_jac(burn_t& state, MatrixType& jac)
{

    // Initialize jac to 0

    jac.zero();

    // Get the data from the state

    const Real rho  = state.rho;
    const Real temp = state.T;
    const Real abar = state.abar;
    const Real zbar = state.zbar;

    Array1D<Real, 1, NumSpec> y;
    for (int i = 1; i <= NumSpec; ++i) {
        y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    const RateEvaluator rr(rho, temp, y);

    // Species Jacobian elements with respect to other species
    // and with respect to temperature

    dFydy_dFydT_aprox13(y, rr, jac);

    // Energy generate rate Jacobian element with respect to temperature

    auto jac_species_temp = [&](int i) -> Real& { 
        return jac(i, net_itemp);
    };

    jac(net_ienuc,net_itemp) = ener_gener_rate(jac_species_temp);

    // Account for the thermal neutrino losses

    Real b1, sneut, dsneutdt, dsneutdd, snuda, snudz;

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    for (int j = 1; j <= NumSpec; ++j) {
       b1 = (-abar * abar * snuda + (zion[j-1] - zbar) * abar * snudz);
       jac.add(net_ienuc, j, -b1);
    }

    jac(net_ienuc,net_itemp) -= dsneutdt;

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc,j));
    }

    // Temperature Jacobian elements

    temperature_jac(state, jac);

}


AMREX_INLINE
void set_up_screening_factors()
{
    // Compute and store the more expensive screening factors

    using namespace Species;

    int jscr = 0;

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], 4.0e0_rt, 8.0e0_rt);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[Ne20-1], aion[Ne20-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[Mg24-1], aion[Mg24-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 13.0e0_rt, 27.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Si28-1], aion[Si28-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 15.0e0_rt, 31.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[S32-1], aion[S32-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 17.0e0_rt, 35.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Ar36-1], aion[Ar36-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 19.0e0_rt, 39.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Ca40-1], aion[Ca40-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 21.0e0_rt, 43.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Ti44-1], aion[Ti44-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 23.0e0_rt, 47.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Cr48-1], aion[Cr48-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 25.0e0_rt, 51.0e0_rt, 1.0e0_rt, 1.0e0_rt);

    add_screening_factor(jscr++, zion[Fe52-1], aion[Fe52-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, 27.0e0_rt, 55.0e0_rt, 1.0e0_rt, 1.0e0_rt);

}

#endif
