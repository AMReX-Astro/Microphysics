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


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dfdy_isotopes_aprox13(Array1D<Real, 1, NumSpec> const& y,
                           burn_t const& state, RateEvaluator& rr,
                           MatrixType& jac)
{

    // this routine sets up the dense aprox13 jacobian for the isotopes

    using namespace Species;
    using namespace Rates;

    // he4 jacobian elements
    // d(he4)/d(he4)
    Array1D<Real, 1, 20> b20 {
        -1.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a),
        -y(C12)  * rr.rates(1,ircag),
        -y(O16)  * rr.rates(1,iroag),
        -y(Ne20) * rr.rates(1,irneag),
        -y(Mg24) * rr.rates(1,irmgag),
        -y(Si28) * rr.rates(1,irsiag),
        -y(S32)  * rr.rates(1,irsag),
        -y(Ar36) * rr.rates(1,irarag),
        -y(Ca40) * rr.rates(1,ircaag),
        -y(Ti44) * rr.rates(1,irtiag),
        -y(Cr48) * rr.rates(1,ircrag),
        -y(Fe52) * rr.rates(1,irfeag),
        -y(Mg24) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1)),
        -y(Si28) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1)),
        -y(S32)  * rr.rates(1,irsap)  * (1.0e0_rt-rr.rates(1,irt1)),
        -y(Ar36) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1)),
        -y(Ca40) * rr.rates(1,ircaap) * (1.0e0_rt-rr.rates(1,irv1)),
        -y(Ti44) * rr.rates(1,irtiap) * (1.0e0_rt-rr.rates(1,irw1)),
        -y(Cr48) * rr.rates(1,ircrap) * (1.0e0_rt-rr.rates(1,irx1)),
        -y(Fe52) * rr.rates(1,irfeap) * (1.0e0_rt-rr.rates(1,iry1))
    };

    jac(He4,He4) = esum20(b20);

    Array1D<Real, 1, 5> b;

    // (he4)/d(o16)
    b(1) =  0.5e0_rt * y(C12) * rr.rates(1,ir1216),
    b(2) =  1.12e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,ir1616),
    b(3) =  0.68e0_rt * rr.rates(1,irs1) * 0.5e0_rt*y(O16) * rr.rates(1,ir1616),
    b(4) =  rr.rates(1,iroga),
    b(5) = -y(He4) * rr.rates(1,iroag),

    jac(He4,O16) = esum5(b);

    // d(he4)/d(c12)
    b(1) =  y(C12) * rr.rates(1,ir1212);
    b(2) =  0.5e0_rt * y(O16) * rr.rates(1,ir1216);
    b(3) =  3.0e0_rt * rr.rates(1,irg3a);
    b(4) = -y(He4) * rr.rates(1,ircag);

    jac(He4,C12) = esum4(b);

    // d(he4)/d(si28)
    b(1) =  rr.rates(1,irsiga);
    b(2) = -y(He4) * rr.rates(1,irsiag);
    b(3) = -y(He4) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));
    b(4) =  rr.rates(1,irr1) * rr.rates(1,irsigp);

    jac(He4,Si28) = esum4(b);

    // d(he4)/d(s32)
    b(1) =  rr.rates(1,irsga);
    b(2) = -y(He4) * rr.rates(1,irsag);
    b(3) = -y(He4) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));
    b(4) =  rr.rates(1,irs1) * rr.rates(1,irsgp);

    jac(He4,S32) = esum4(b);

    // d(he4)/d(ar36)
    b(1)  =  rr.rates(1,irarga);
    b(2)  = -y(He4) * rr.rates(1,irarag);
    b(3)  = -y(He4) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1));
    b(4)  =  rr.rates(1,irt1) * rr.rates(1,irargp);

    jac(He4,Ar36) = esum4(b);

    // d(he4)/d(ca40)
    b(1) =  rr.rates(1,ircaga);
    b(2) = -y(He4) * rr.rates(1,ircaag);
    b(3) = -y(He4) * rr.rates(1,ircaap) * (1.0e0_rt-rr.rates(1,irv1));
    b(4) =  rr.rates(1,iru1) * rr.rates(1,ircagp);

    jac(He4,Ca40) = esum4(b);

    // d(he4)/d(ti44)
    b(1) =  rr.rates(1,irtiga);
    b(2) = -y(He4) * rr.rates(1,irtiag);
    b(3) = -y(He4) * rr.rates(1,irtiap) * (1.0e0_rt-rr.rates(1,irw1));
    b(4) =  rr.rates(1,irv1) * rr.rates(1,irtigp);

    jac(He4,Ti44) = esum4(b);

    // d(he4)/d(cr48)
    b(1) =  rr.rates(1,ircrga);
    b(2) = -y(He4) * rr.rates(1,ircrag);
    b(3) = -y(He4) * rr.rates(1,ircrap) * (1.0e0_rt-rr.rates(1,irx1));
    b(4) =  rr.rates(1,irw1) * rr.rates(1,ircrgp);

    jac(He4,Cr48) = esum4(b);

    // d(he4)/d(fe52)
    b(1) =  rr.rates(1,irfega);
    b(2) = -y(He4) * rr.rates(1,irfeag);
    b(3) = -y(He4) * rr.rates(1,irfeap) * (1.0e0_rt-rr.rates(1,iry1));
    b(4) =  rr.rates(1,irx1) * rr.rates(1,irfegp);

    jac(He4,Fe52) = esum4(b);

    // d(c12)/d(c12)
    b(1) = -2.0e0_rt * y(C12) * rr.rates(1,ir1212);
    b(2) = -y(O16) * rr.rates(1,ir1216);
    b(3) = -rr.rates(1,irg3a);
    b(4) = -y(He4) * rr.rates(1,ircag);

    jac(C12,C12) = esum4(b);

    // d(o16)/d(o16)
    b(1) = -y(C12) * rr.rates(1,ir1216);
    b(2) = -2.0e0_rt * y(O16) * rr.rates(1,ir1616);
    b(3) = -y(He4) * rr.rates(1,iroag);
    b(4) = -rr.rates(1,iroga);

    jac(O16,O16) = esum4(b);

    // si28 jacobian elements
    // d(si28)/d(he4)
    b(1) =  y(Mg24) * rr.rates(1,irmgag);
    b(2) = -y(Si28) * rr.rates(1,irsiag);
    b(3) =  y(Mg24) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));
    b(4) = -y(Si28) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));

    jac(Si28,He4) = esum4(b);

    // d(si28)/d(si28)
    b(1) =  -y(He4) * rr.rates(1,irsiag);
    b(2) = -rr.rates(1,irsiga);
    b(3) = -rr.rates(1,irr1) * rr.rates(1,irsigp);
    b(4) = -y(He4) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));

    jac(Si28,Si28) = esum4(b);

    // s32 jacobian elements
    // d(s32)/d(he4)
    b(1) =  y(Si28) * rr.rates(1,irsiag);
    b(2) = -y(S32) * rr.rates(1,irsag);
    b(3) =  y(Si28) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));
    b(4) = -y(S32) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));

    jac(S32,He4) = esum4(b);

    // d(s32)/d(s32)
    b(1) = -y(He4) * rr.rates(1,irsag);
    b(2) = -rr.rates(1,irsga);
    b(3) = -rr.rates(1,irs1) * rr.rates(1,irsgp);
    b(4) = -y(He4) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));

    jac(S32,S32) = esum4(b);

    // ar36 jacobian elements
    // d(ar36)/d(he4)
    b(1) =  y(S32)  * rr.rates(1,irsag);
    b(2) = -y(Ar36) * rr.rates(1,irarag);
    b(3) =  y(S32)  * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));
    b(4) = -y(Ar36) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1));

    jac(Ar36,He4) = esum4(b);

    // d(ar36)/d(ar36)
    b(1) = -y(He4) * rr.rates(1,irarag);
    b(2) = -rr.rates(1,irarga);
    b(3) = -rr.rates(1,irt1) * rr.rates(1,irargp);
    b(4) = -y(He4) * rr.rates(1,irarap) * (1.0e0_rt-rr.rates(1,iru1));

    jac(Ar36,Ar36) = esum4(b);

    // ca40 jacobian elements
    // d(ca40)/d(he4)
    b(1)  =  y(Ar36) * rr.rates(1,irarag);
    b(2)  = -y(Ca40) * rr.rates(1,ircaag);
    b(3)  =  y(Ar36) * rr.rates(1,irarap)*(1.0e0_rt-rr.rates(1,iru1));
    b(4)  = -y(Ca40) * rr.rates(1,ircaap)*(1.0e0_rt-rr.rates(1,irv1));

    jac(Ca40,He4) = esum4(b);

    // d(ca40)/d(ca40)
    b(1) = -y(He4) * rr.rates(1, ircaag);
    b(2) = -rr.rates(1, ircaga);
    b(3) = -rr.rates(1, ircagp) * rr.rates(1, iru1);
    b(4) = -y(He4) * rr.rates(1, ircaap)*(1.0e0_rt-rr.rates(1, irv1));

    jac(Ca40,Ca40) = esum4(b);

    // ti44 jacobian elements
    // d(ti44)/d(he4)
    b(1) =  y(Ca40) * rr.rates(1, ircaag);
    b(2) = -y(Ti44) * rr.rates(1, irtiag);
    b(3) =  y(Ca40) * rr.rates(1, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
    b(4) = -y(Ti44) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));

    jac(Ti44,He4) = esum4(b);

    // d(ti44)/d(ti44)
    b(1) = -y(He4) * rr.rates(1, irtiag);
    b(2) = -rr.rates(1, irtiga);
    b(3) = -rr.rates(1, irv1) * rr.rates(1, irtigp);
    b(4) = -y(He4) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));

    jac(Ti44,Ti44) = esum4(b);

    // cr48 jacobian elements
    // d(cr48)/d(he4)
    b(1) =  y(Ti44) * rr.rates(1, irtiag);
    b(2) = -y(Cr48) * rr.rates(1, ircrag);
    b(3) =  y(Ti44) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));
    b(4) = -y(Cr48) * rr.rates(1, ircrap)*(1.0e0_rt-rr.rates(1, irx1));

    jac(Cr48,He4) = esum4(b);

    // d(cr48)/d(cr48)
    b(1) = -y(He4) * rr.rates(1, ircrag);
    b(2) = -rr.rates(1, ircrga);
    b(3) = -rr.rates(1, irw1) * rr.rates(1, ircrgp);
    b(4) = -y(He4) * rr.rates(1, ircrap)*(1.0e0_rt-rr.rates(1, irx1));

    jac(Cr48,Cr48) = esum4(b);

    // fe52 jacobian elements
    // d(fe52)/d(he4)
    b(1) =  y(Cr48) * rr.rates(1, ircrag);
    b(2) = -y(Fe52) * rr.rates(1, irfeag);
    b(3) =  y(Cr48) * rr.rates(1, ircrap) * (1.0e0_rt-rr.rates(1, irx1));
    b(4) = -y(Fe52) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

    jac(Fe52,He4) = esum4(b);

    // d(fe52)/d(fe52)
    b(1) = -y(He4) * rr.rates(1, irfeag);
    b(2) = -rr.rates(1, irfega);
    b(3) = -rr.rates(1, irx1) * rr.rates(1, irfegp);
    b(4) = -y(He4) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

    jac(Fe52,Fe52) = esum4(b);

    // d(he4)/d(mg24)
    b(1) =  rr.rates(1,irmgga);
    b(2) = -y(He4) * rr.rates(1,irmgag);
    b(3) = -y(He4) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

    jac(He4,Mg24) = esum3(b);

    // mg24 jacobian elements
    // d(mg24)/d(he4)
    b(1) =  y(Ne20) * rr.rates(1,irneag);
    b(2) = -y(Mg24) * rr.rates(1,irmgag);
    b(3) = -y(Mg24) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

    jac(Mg24,He4) = esum3(b);

    // d(mg24)/d(mg24)
    b(1) = -y(He4) * rr.rates(1,irmgag);
    b(2) = -rr.rates(1,irmgga);
    b(3) = -y(He4) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

    jac(Mg24,Mg24) = esum3(b);

    // d(si28)/d(o16)
    b(1) = 0.5e0_rt * y(C12) * rr.rates(1,ir1216);
    b(2) = 1.12e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,ir1616);
    b(3) = 0.68e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,irs1) * rr.rates(1,ir1616);

    jac(Si28,O16) = esum3(b);

    // d(he4)/d(ne20)
    b(1) =  rr.rates(1,irnega);
    b(2) = -y(He4) * rr.rates(1,irneag);

    jac(He4,Ne20) = b(1) + b(2);

    // d(he4)/d(ni56)
    b(1) = rr.rates(1,irniga);
    b(2) = rr.rates(1,iry1) * rr.rates(1,irnigp);

    jac(He4,Ni56) = b(1) + b(2);

    // c12 jacobian elements
    // d(c12)/d(he4)
    b(1) =  0.5e0_rt * y(He4) * y(He4) * rr.rates(1,ir3a);
    b(2) = -y(C12) * rr.rates(1,ircag);

    jac(C12,He4) = b(1) + b(2);

    // d(c12)/d(o16)
    b(1) = -y(C12) * rr.rates(1,ir1216);
    b(2) =  rr.rates(1,iroga);

    jac(C12,O16) = b(1) + b(2);

    // o16 jacobian elements
    // d(o16)/d(he4)
    b(1) =  y(C12)*rr.rates(1,ircag);
    b(2) = -y(O16)*rr.rates(1,iroag);

    jac(O16,He4) = b(1) + b(2);

    // d(o16)/d(c12)
    b(1) = -y(O16)*rr.rates(1,ir1216);
    b(2) =  y(He4)*rr.rates(1,ircag);

    jac(O16,C12) = b(1) + b(2);

    // ne20 jacobian elements
    // d(ne20)/d(he4)
    b(1) =  y(O16) * rr.rates(1,iroag);
    b(2) = -y(Ne20) * rr.rates(1,irneag);

    jac(Ne20,He4) = b(1) + b(2);

    // d(ne20)/d(ne20)
    b(1) = -y(He4) * rr.rates(1,irneag);
    b(2) = -rr.rates(1,irnega);

    jac(Ne20,Ne20) = b(1) + b(2);

    // d(mg24)/d(si28)
    b(1) = rr.rates(1,irsiga);
    b(2) = rr.rates(1,irr1) * rr.rates(1,irsigp);

    jac(Mg24,Si28) = b(1) + b(2);

    // d(si28)/d(mg24)
    b(1) =  y(He4) * rr.rates(1,irmgag);
    b(2) =  y(He4) * rr.rates(1,irmgap) * (1.0e0_rt-rr.rates(1,irr1));

    jac(Si28,Mg24) = b(1) + b(2);

    // d(si28)/d(s32)
    b(1) = rr.rates(1,irsga);
    b(2) = rr.rates(1,irs1) * rr.rates(1,irsgp);

    jac(Si28,S32) = b(1) + b(2);

    // d(s32)/d(o16)
    b(1) = 0.68e0_rt*0.5e0_rt*y(O16)*rr.rates(1,ir1616)*(1.0e0_rt-rr.rates(1,irs1));
    b(2) = 0.2e0_rt * 0.5e0_rt*y(O16) * rr.rates(1,ir1616);

    jac(S32,O16) = b(1) + b(2);

    // d(s32)/d(si28)
    b(1)  =y(He4) * rr.rates(1,irsiag);
    b(2) = y(He4) * rr.rates(1,irsiap) * (1.0e0_rt-rr.rates(1,irs1));

    jac(S32,Si28) = b(1) + b(2);

    // d(s32)/d(ar36)
    b(1) = rr.rates(1,irarga);
    b(2) = rr.rates(1,irt1) * rr.rates(1,irargp);

    jac(S32,Ar36) = b(1) + b(2);

    // d(ar36)/d(s32)
    b(1) = y(He4) * rr.rates(1,irsag);
    b(2) = y(He4) * rr.rates(1,irsap) * (1.0e0_rt-rr.rates(1,irt1));

    jac(Ar36,S32) = b(1) + b(2);

    // d(ar36)/d(ca40)
    b(1) = rr.rates(1,ircaga);
    b(2) = rr.rates(1,ircagp) * rr.rates(1,iru1);

    jac(Ar36,Ca40) = b(1) + b(2);

    // d(ca40)/d(ar36)
    b(1) =  y(He4) * rr.rates(1,irarag);
    b(2) =  y(He4) * rr.rates(1,irarap)*(1.0e0_rt-rr.rates(1,iru1));

    jac(Ca40,Ar36) = b(1) + b(2);

    // d(ca40)/d(ti44)
    b(1) = rr.rates(1, irtiga);
    b(2) = rr.rates(1, irtigp) * rr.rates(1, irv1);

    jac(Ca40,Ti44) = b(1) + b(2);

    // d(ti44)/d(ca40)
    b(1) =  y(He4) * rr.rates(1, ircaag);
    b(2) =  y(He4) * rr.rates(1, ircaap)*(1.0e0_rt-rr.rates(1, irv1));

    jac(Ti44,Ca40) = b(1) + b(2);

    // d(ti44)/d(cr48)
    b(1) = rr.rates(1, ircrga);
    b(2) = rr.rates(1, irw1) * rr.rates(1, ircrgp);

    jac(Ti44,Cr48) = b(1) + b(2);

    // d(cr48)/d(ti44)
    b(1) =  y(He4) * rr.rates(1, irtiag);
    b(2) =  y(He4) * rr.rates(1, irtiap)*(1.0e0_rt-rr.rates(1, irw1));

    jac(Cr48,Ti44) = b(1) + b(2);

    // d(cr48)/d(fe52)
    b(1) = rr.rates(1, irfega);
    b(2) = rr.rates(1, irx1) * rr.rates(1, irfegp);

    jac(Cr48,Fe52) = b(1) + b(2);

    // d(fe52)/d(cr48)
    b(1) = y(He4) * rr.rates(1, ircrag);
    b(2) = y(He4) * rr.rates(1, ircrap) * (1.0e0_rt-rr.rates(1, irx1));

    jac(Fe52,Cr48) = b(1) + b(2);

    // d(fe52)/d(ni56)
    b(1) = rr.rates(1, irniga);
    b(2) = rr.rates(1, iry1) * rr.rates(1, irnigp);

    jac(Fe52,Ni56) = b(1) + b(2);

    // ni56 jacobian elements
    // d(ni56)/d(he4)
    b(1) =  y(Fe52) * rr.rates(1, irfeag);
    b(2) =  y(Fe52) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

    jac(Ni56,He4) = b(1) + b(2);

    // d(ni56)/d(fe52)
    b(1) = y(He4) * rr.rates(1, irfeag);
    b(2) = y(He4) * rr.rates(1, irfeap) * (1.0e0_rt-rr.rates(1, iry1));

    jac(Ni56,Fe52) = b(1) + b(2);

    // d(ni56)/d(ni56)
    b(1) = -rr.rates(1, irniga);
    b(2) = -rr.rates(1, iry1) * rr.rates(1, irnigp);

        jac(Ni56,Ni56) = b(1) + b(2);

    // d(o16)/d(ne20)
    jac(O16,Ne20) = rr.rates(1,irnega);

    // d(ne20)/d(c12)
    jac(Ne20,C12) = y(C12) * rr.rates(1,ir1212);

    // d(ne20)/d(o16)
    jac(Ne20,O16) = y(He4) * rr.rates(1,iroag);

    // d(ne20)/d(mg24)
    jac(Ne20,Mg24) = rr.rates(1,irmgga);

    // d(mg24)/d(c12)
    jac(Mg24,C12) = 0.5e0_rt * y(O16) * rr.rates(1,ir1216);

    // d(mg24)/d(o16)
    jac(Mg24,O16) = 0.5e0_rt * y(C12) * rr.rates(1,ir1216);

    // d(mg24)/d(ne20)
    jac(Mg24,Ne20) = y(He4) * rr.rates(1,irneag);

    // d(si28)/d(c12)
    jac(Si28,C12) = 0.5e0_rt * y(O16) * rr.rates(1,ir1216);

}

template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, Real& enuc)
{

    using namespace aprox13;

    // Computes the instantaneous energy generation rate

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives
    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2
    enuc = Xdot * C::Legacy::enuc_conv2;

}


// Evaluates the right hand side of the aprox13 ODEs
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs(Array1D<Real, 1, NumSpec> const& y, RateEvaluator const& rr,
         Array1D<Real, 1, neqs>& dydt,
         const bool deriva, const bool for_jacobian_tderiv)
{
    using namespace Species;
    using namespace Rates;

    // deriva is used in forming the analytic Jacobian to get
    // the derivative wrt A

    const int index_rate = for_jacobian_tderiv ? 2 : 1;

    for (int i = 1; i <= NumSpec; ++i) {
        dydt(i) = 0.0_rt;
    }

    Array1D<Real, 1, 17> a;

    // he4 reactions
    // heavy ion reactions
    a(1)  = 0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate, ir1212);
    a(2)  = 0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(3)  = 0.56e0_rt * 0.5e0_rt * y(O16) * y(O16) * rr.rates(index_rate, ir1616);

    dydt(He4) = dydt(He4) + esum3(a);

    // (a,g) and (g,a) reactions
    a(1)  = -0.5e0_rt * y(He4) * y(He4) * y(He4) * rr.rates(index_rate, ir3a);
    a(2)  =  3.0e0_rt * y(C12) * rr.rates(index_rate, irg3a);
    a(3)  = -y(He4)  * y(C12) * rr.rates(index_rate, ircag);
    a(4)  =  y(O16)  * rr.rates(index_rate, iroga);
    a(5)  = -y(He4)  * y(O16) * rr.rates(index_rate, iroag);
    a(6)  =  y(Ne20) * rr.rates(index_rate, irnega);
    a(7)  = -y(He4)  * y(Ne20) * rr.rates(index_rate, irneag);
    a(8)  =  y(Mg24) * rr.rates(index_rate, irmgga);
    a(9)  = -y(He4)  * y(Mg24)* rr.rates(index_rate, irmgag);
    a(10) =  y(Si28) * rr.rates(index_rate, irsiga);
    a(11) = -y(He4)  * y(Si28)*rr.rates(index_rate, irsiag);
    a(12) =  y(S32)  * rr.rates(index_rate, irsga);

    dydt(He4) = dydt(He4) + esum12(a);

    a(1)  = -y(He4)  * y(S32) * rr.rates(index_rate, irsag);
    a(2)  =  y(Ar36) * rr.rates(index_rate, irarga);
    a(3)  = -y(He4)  * y(Ar36)*rr.rates(index_rate, irarag);
    a(4)  =  y(Ca40) * rr.rates(index_rate, ircaga);
    a(5)  = -y(He4)  * y(Ca40)*rr.rates(index_rate, ircaag);
    a(6)  =  y(Ti44) * rr.rates(index_rate, irtiga);
    a(7)  = -y(He4)  * y(Ti44)*rr.rates(index_rate, irtiag);
    a(8)  =  y(Cr48) * rr.rates(index_rate, ircrga);
    a(9)  = -y(He4)  * y(Cr48)*rr.rates(index_rate, ircrag);
    a(10) =  y(Fe52) * rr.rates(index_rate, irfega);
    a(11) = -y(He4)  * y(Fe52) * rr.rates(index_rate, irfeag);
    a(12) =  y(Ni56) * rr.rates(index_rate, irniga);

    dydt(He4) = dydt(He4) + esum12(a);

    // (a,p)(p,g) and (g,p)(p,a) reactions

    if (!deriva) {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*rr.rates(index_rate, irs1)*rr.rates(index_rate, ir1616);
       a(2)  = -y(He4)  * y(Mg24) * rr.rates(index_rate, irmgap)*(1.0e0_rt-rr.rates(index_rate, irr1));
       a(3)  =  y(Si28) * rr.rates(index_rate, irsigp) * rr.rates(index_rate, irr1);
       a(4)  = -y(He4)  * y(Si28) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(5)  =  y(S32)  * rr.rates(index_rate, irsgp) * rr.rates(index_rate, irs1);
       a(6)  = -y(He4)  * y(S32) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(index_rate, irt1));
       a(7)  =  y(Ar36) * rr.rates(index_rate, irargp) * rr.rates(index_rate, irt1);
       a(8)  = -y(He4)  * y(Ar36) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(index_rate, iru1));
       a(9)  =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(index_rate, iru1);
       a(10) = -y(He4)  * y(Ca40) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(index_rate, irv1));
       a(11) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(index_rate, irv1);
       a(12) = -y(He4)  * y(Ti44) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(index_rate, irw1));
       a(13) =  y(Cr48) * rr.rates(index_rate, ircrgp) * rr.rates(index_rate, irw1);
       a(14) = -y(He4)  * y(Cr48) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(index_rate, irx1));
       a(15) =  y(Fe52) * rr.rates(index_rate, irfegp) * rr.rates(index_rate, irx1);
       a(16) = -y(He4)  * y(Fe52) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(index_rate, iry1));
       a(17) =  y(Ni56) * rr.rates(index_rate, irnigp) * rr.rates(index_rate, iry1);

       dydt(He4) = dydt(He4) + esum17(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(1, irs1) * rr.rates(index_rate, ir1616);
       a(2)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(index_rate, irs1) * rr.rates(1, ir1616);
       a(3)  = -y(He4)*y(Mg24) * rr.rates(index_rate, irmgap)*(1.0e0_rt - rr.rates(1, irr1));
       a(4)  =  y(He4)*y(Mg24) * rr.rates(1, irmgap)*rr.rates(index_rate, irr1);
       a(5)  =  y(Si28) * rr.rates(1, irsigp) * rr.rates(index_rate, irr1);
       a(6)  =  y(Si28) * rr.rates(index_rate, irsigp) * rr.rates(1, irr1);
       a(7)  = -y(He4)*y(Si28) * rr.rates(index_rate, irsiap)*(1.0e0_rt - rr.rates(1, irs1));
       a(8)  =  y(He4)*y(Si28) * rr.rates(1, irsiap) * rr.rates(index_rate, irs1);
       a(9)  =  y(S32)  * rr.rates(1, irsgp) * rr.rates(index_rate, irs1);
       a(10) =  y(S32)  * rr.rates(index_rate, irsgp) * rr.rates(1, irs1);

       dydt(He4) = dydt(He4) + esum10(a);

       a(1)  = -y(He4)*y(S32) * rr.rates(index_rate, irsap)*(1.0e0_rt - rr.rates(1, irt1));
       a(2)  =  y(He4)*y(S32) * rr.rates(1, irsap)*rr.rates(index_rate, irt1);
       a(3)  =  y(Ar36) * rr.rates(1, irargp) * rr.rates(index_rate, irt1);
       a(4)  =  y(Ar36) * rr.rates(index_rate, irargp) * rr.rates(1, irt1);
       a(5)  = -y(He4)*y(Ar36) * rr.rates(index_rate, irarap)*(1.0e0_rt - rr.rates(1, iru1));
       a(6)  =  y(He4)*y(Ar36) * rr.rates(1, irarap)*rr.rates(index_rate, iru1);
       a(7)  =  y(Ca40) * rr.rates(1, ircagp) * rr.rates(index_rate, iru1);
       a(8)  =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(1, iru1);
       a(9)  = -y(He4)*y(Ca40) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
       a(10) =  y(He4)*y(Ca40) * rr.rates(1, ircaap)*rr.rates(index_rate, irv1);
       a(11) =  y(Ti44) * rr.rates(1, irtigp) * rr.rates(index_rate, irv1);
       a(12) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(1, irv1);

       dydt(He4) = dydt(He4) + esum12(a);

       a(1)  = -y(He4)*y(Ti44) * rr.rates(index_rate, irtiap)*(1.0e0_rt - rr.rates(1, irw1));
       a(2)  =  y(He4)*y(Ti44) * rr.rates(1, irtiap)*rr.rates(index_rate, irw1);
       a(3)  =  y(Cr48) * rr.rates(1, ircrgp) * rr.rates(index_rate, irw1);
       a(4)  =  y(Cr48) * rr.rates(index_rate, ircrgp) * rr.rates(1, irw1);
       a(5)  = -y(He4)*y(Cr48) * rr.rates(index_rate, ircrap)*(1.0e0_rt - rr.rates(1, irx1));
       a(6)  =  y(He4)*y(Cr48) * rr.rates(1, ircrap)*rr.rates(index_rate, irx1);
       a(7)  =  y(Fe52) * rr.rates(1, irfegp) * rr.rates(index_rate, irx1);
       a(8)  =  y(Fe52) * rr.rates(index_rate, irfegp) * rr.rates(1, irx1);
       a(9)  = -y(He4)*y(Fe52) * rr.rates(index_rate, irfeap)*(1.0e0_rt - rr.rates(1, iry1));
       a(10) =  y(He4)*y(Fe52) * rr.rates(1, irfeap)*rr.rates(index_rate, iry1);
       a(11) =  y(Ni56) * rr.rates(1, irnigp) * rr.rates(index_rate, iry1);
       a(12) =  y(Ni56) * rr.rates(index_rate, irnigp) * rr.rates(1, iry1);

       dydt(He4) = dydt(He4) + esum12(a);
    
    }

    // c12 reactions
    a(1) = -y(C12) * y(C12) * rr.rates(index_rate, ir1212);
    a(2) = -y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(3) =  y(He4) * y(He4) * y(He4) * rr.rates(index_rate, ir3a) / 6.0_rt;;
    a(4) = -y(C12) * rr.rates(index_rate, irg3a);
    a(5) = -y(C12) * y(He4) * rr.rates(index_rate, ircag);
    a(6) =  y(O16) * rr.rates(index_rate, iroga);

    dydt(C12) = dydt(C12) + esum6(a);


    // o16 reactions
    a(1) = -y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(2) = -y(O16) * y(O16) * rr.rates(index_rate, ir1616);
    a(3) =  y(C12) * y(He4) * rr.rates(index_rate, ircag);
    a(4) = -y(O16) * y(He4) * rr.rates(index_rate, iroag);
    a(5) = -y(O16) * rr.rates(index_rate, iroga);
    a(6) =  y(Ne20) * rr.rates(index_rate, irnega);

    dydt(O16) = dydt(O16) + esum6(a);


    // ne20 reactions
    a(1) =  0.5e0_rt * y(C12) * y(C12) * rr.rates(index_rate, ir1212);
    a(2) =  y(O16) * y(He4) * rr.rates(index_rate, iroag);
    a(3) = -y(Ne20) * y(He4) * rr.rates(index_rate, irneag);
    a(4) = -y(Ne20) * rr.rates(index_rate, irnega);
    a(5) =  y(Mg24) * rr.rates(index_rate, irmgga);

    dydt(Ne20) = dydt(Ne20) + esum5(a);


    // mg24 reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(2) =  y(Ne20) * y(He4) * rr.rates(index_rate, irneag);
    a(3) = -y(Mg24) * y(He4) * rr.rates(index_rate, irmgag);
    a(4) = -y(Mg24) * rr.rates(index_rate, irmgga);
    a(5) =  y(Si28) * rr.rates(index_rate, irsiga);

    dydt(Mg24) = dydt(Mg24) + esum5(a);

    if (!deriva) {

       a(1) = -y(Mg24) * y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt-rr.rates(index_rate, irr1));
       a(2) =  y(Si28) * rr.rates(index_rate, irr1) * rr.rates(index_rate, irsigp);

       dydt(Mg24) = dydt(Mg24) + a(1) + a(2);

    } else {

       a(1) = -y(Mg24)*y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt - rr.rates(1, irr1));
       a(2) =  y(Mg24)*y(He4) * rr.rates(1, irmgap)*rr.rates(index_rate, irr1);
       a(3) =  y(Si28) * rr.rates(1, irr1) * rr.rates(index_rate, irsigp);
       a(4) =  y(Si28) * rr.rates(index_rate, irr1) * rr.rates(1, irsigp);

       dydt(Mg24) = dydt(Mg24) + esum4(a);

    }


    // si28 reactions
    a(1) =  0.5e0_rt * y(C12) * y(O16) * rr.rates(index_rate, ir1216);
    a(2) =  0.56e0_rt * 0.5e0_rt*y(O16) * y(O16) * rr.rates(index_rate, ir1616);
    a(3) =  y(Mg24) * y(He4) * rr.rates(index_rate, irmgag);
    a(4) = -y(Si28) * y(He4) * rr.rates(index_rate, irsiag);
    a(5) = -y(Si28) * rr.rates(index_rate, irsiga);
    a(6) =  y(S32)  * rr.rates(index_rate, irsga);

    dydt(Si28) = dydt(Si28) + esum6(a);

    if (!deriva) {

       a(1) =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)*rr.rates(index_rate, irs1)*rr.rates(index_rate, ir1616);
       a(2) =  y(Mg24) * y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt-rr.rates(index_rate, irr1));
       a(3) = -y(Si28) * rr.rates(index_rate, irr1) * rr.rates(index_rate, irsigp);
       a(4) = -y(Si28) * y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(5) =  y(S32)  * rr.rates(index_rate, irs1) * rr.rates(index_rate, irsgp);

       dydt(Si28) = dydt(Si28) + esum5(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(1, irs1)*rr.rates(index_rate, ir1616);
       a(2)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(index_rate, irs1)*rr.rates(1, ir1616);
       a(3)  =  y(Mg24)*y(He4) * rr.rates(index_rate, irmgap)*(1.0e0_rt - rr.rates(1, irr1));
       a(4)  = -y(Mg24)*y(He4) * rr.rates(1, irmgap)*rr.rates(index_rate, irr1);
       a(5)  = -y(Si28) * rr.rates(1, irr1) * rr.rates(index_rate, irsigp);
       a(6)  = -y(Si28) * rr.rates(index_rate, irr1) * rr.rates(1, irsigp);
       a(7)  = -y(Si28)*y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt - rr.rates(1, irs1));
       a(8)  =  y(Si28)*y(He4) * rr.rates(1, irsiap)*rr.rates(index_rate, irs1);
       a(9)  = y(S32) * rr.rates(1, irs1) * rr.rates(index_rate, irsgp);
       a(10) = y(S32) * rr.rates(index_rate, irs1) * rr.rates(1, irsgp);

       dydt(Si28) = dydt(Si28) + esum10(a);

    }


    // s32 reactions
    a(1) =  0.1e0_rt * 0.5e0_rt*y(O16) * y(O16) * rr.rates(index_rate, ir1616);
    a(2) =  y(Si28) * y(He4) * rr.rates(index_rate, irsiag);
    a(3) = -y(S32) * y(He4) * rr.rates(index_rate, irsag);
    a(4) = -y(S32) * rr.rates(index_rate, irsga);
    a(5) =  y(Ar36) * rr.rates(index_rate, irarga);

    dydt(S32) = dydt(S32) + esum5(a);


    if (!deriva) {

       a(1) =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16)* rr.rates(index_rate, ir1616)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(2) =  y(Si28) * y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(index_rate, irs1));
       a(3) = -y(S32) * rr.rates(index_rate, irs1) * rr.rates(index_rate, irsgp);
       a(4) = -y(S32) * y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(index_rate, irt1));
       a(5) =  y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(index_rate, irargp);

       dydt(S32) = dydt(S32) + esum5(a);

    } else {

       a(1)  =  0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(index_rate, ir1616)*(1.0e0_rt-rr.rates(1, irs1));
       a(2)  = -0.34e0_rt*0.5e0_rt*y(O16)*y(O16) * rr.rates(1, ir1616)*rr.rates(index_rate, irs1);
       a(3)  =  y(Si28)*y(He4) * rr.rates(index_rate, irsiap)*(1.0e0_rt-rr.rates(1, irs1));
       a(4)  = -y(Si28)*y(He4) * rr.rates(1, irsiap)*rr.rates(index_rate, irs1);
       a(5)  = -y(S32) * rr.rates(1, irs1) * rr.rates(index_rate, irsgp);
       a(6)  = -y(S32) * rr.rates(index_rate, irs1) * rr.rates(1, irsgp);
       a(7)  = -y(S32)*y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(1, irt1));
       a(8)  =  y(S32)*y(He4) * rr.rates(1, irsap)*rr.rates(index_rate, irt1);
       a(9)  =  y(Ar36) * rr.rates(1, irt1) * rr.rates(index_rate, irargp);
       a(10) =  y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(1, irargp);

       dydt(S32) = dydt(S32) + esum10(a);

    }


    // ar36 reactions
    a(1) =  y(S32)  * y(He4) * rr.rates(index_rate, irsag);
    a(2) = -y(Ar36) * y(He4) * rr.rates(index_rate, irarag);
    a(3) = -y(Ar36) * rr.rates(index_rate, irarga);
    a(4) =  y(Ca40) * rr.rates(index_rate, ircaga);

    dydt(Ar36) = dydt(Ar36) + esum4(a);

    if (!deriva) {

       a(1) = y(S32)  * y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt-rr.rates(index_rate, irt1));
       a(2) = -y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(index_rate, irargp);
       a(3) = -y(Ar36) * y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(index_rate, iru1));
       a(4) =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(index_rate, iru1);

       dydt(Ar36) = dydt(Ar36) + esum4(a);

    } else {

       a(1) =  y(S32)*y(He4) * rr.rates(index_rate, irsap)*(1.0e0_rt - rr.rates(1, irt1));
       a(2) = -y(S32)*y(He4) * rr.rates(1, irsap)*rr.rates(index_rate, irt1);
       a(3) = -y(Ar36) * rr.rates(1, irt1) * rr.rates(index_rate, irargp);
       a(4) = -y(Ar36) * rr.rates(index_rate, irt1) * rr.rates(1, irargp);
       a(5) = -y(Ar36)*y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(1, iru1));
       a(6) =  y(Ar36)*y(He4) * rr.rates(1, irarap)*rr.rates(index_rate, iru1);
       a(7) =  y(Ca40) * rr.rates(1, ircagp) * rr.rates(index_rate, iru1);
       a(8) =  y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(1, iru1);

       dydt(Ar36) = dydt(Ar36) + esum8(a);

    }


    // ca40 reactions
    a(1) =  y(Ar36) * y(He4) * rr.rates(index_rate, irarag);
    a(2) = -y(Ca40) * y(He4) * rr.rates(index_rate, ircaag);
    a(3) = -y(Ca40) * rr.rates(index_rate, ircaga);
    a(4) =  y(Ti44) * rr.rates(index_rate, irtiga);

    dydt(Ca40) = dydt(Ca40) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ar36) * y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(index_rate, iru1));
       a(2) = -y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(index_rate, iru1);
       a(3) = -y(Ca40) * y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(index_rate, irv1));
       a(4) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(index_rate, irv1);

       dydt(Ca40) = dydt(Ca40) + esum4(a);

    } else {

       a(1) =  y(Ar36)*y(He4) * rr.rates(index_rate, irarap)*(1.0e0_rt-rr.rates(1, iru1));
       a(2) = -y(Ar36)*y(He4) * rr.rates(1, irarap)*rr.rates(index_rate, iru1);
       a(3) = -y(Ca40) * rr.rates(1, ircagp) * rr.rates(index_rate, iru1);
       a(4) = -y(Ca40) * rr.rates(index_rate, ircagp) * rr.rates(1, iru1);
       a(5) = -y(Ca40)*y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
       a(6) =  y(Ca40)*y(He4) * rr.rates(1, ircaap)*rr.rates(index_rate, irv1);
       a(7) =  y(Ti44) * rr.rates(1, irtigp) * rr.rates(index_rate, irv1);
       a(8) =  y(Ti44) * rr.rates(index_rate, irtigp) * rr.rates(1, irv1);

       dydt(Ca40) = dydt(Ca40) + esum8(a);

    }


    // ti44 reactions
    a(1) =  y(Ca40) * y(He4) * rr.rates(index_rate, ircaag);
    a(2) = -y(Ti44) * y(He4) * rr.rates(index_rate, irtiag);
    a(3) = -y(Ti44) * rr.rates(index_rate, irtiga);
    a(4) =  y(Cr48) * rr.rates(index_rate, ircrga);

    dydt(Ti44) = dydt(Ti44) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ca40) * y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(index_rate, irv1));
       a(2) = -y(Ti44) * rr.rates(index_rate, irv1) * rr.rates(index_rate, irtigp);
       a(3) = -y(Ti44) * y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(index_rate, irw1));
       a(4) =  y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(index_rate, ircrgp);

       dydt(Ti44) = dydt(Ti44) + esum4(a);

    } else {

       a(1) =  y(Ca40)*y(He4) * rr.rates(index_rate, ircaap)*(1.0e0_rt-rr.rates(1, irv1));
       a(2) = -y(Ca40)*y(He4) * rr.rates(1, ircaap)*rr.rates(index_rate, irv1);
       a(3) = -y(Ti44) * rr.rates(1, irv1) * rr.rates(index_rate, irtigp);
       a(4) = -y(Ti44) * rr.rates(index_rate, irv1) * rr.rates(1, irtigp);
       a(5) = -y(Ti44)*y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(1, irw1));
       a(6) =  y(Ti44)*y(He4) * rr.rates(1, irtiap)*rr.rates(index_rate, irw1);
       a(7) =  y(Cr48) * rr.rates(1, irw1) * rr.rates(index_rate, ircrgp);
       a(8) =  y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(1, ircrgp);

       dydt(Ti44) = dydt(Ti44) + esum8(a);

    }


    // cr48 reactions
    a(1) =  y(Ti44) * y(He4) * rr.rates(index_rate, irtiag);
    a(2) = -y(Cr48) * y(He4) * rr.rates(index_rate, ircrag);
    a(3) = -y(Cr48) * rr.rates(index_rate, ircrga);
    a(4) =  y(Fe52) * rr.rates(index_rate, irfega);

    dydt(Cr48) = dydt(Cr48) + esum4(a);

    if (!deriva) {

       a(1) =  y(Ti44) * y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(index_rate, irw1));
       a(2) = -y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(index_rate, ircrgp);
       a(3) = -y(Cr48) * y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(index_rate, irx1));
       a(4) =  y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(index_rate, irfegp);

       dydt(Cr48) = dydt(Cr48) + esum4(a);

    } else {

       a(1) =  y(Ti44)*y(He4) * rr.rates(index_rate, irtiap)*(1.0e0_rt-rr.rates(1, irw1));
       a(2) = -y(Ti44)*y(He4) * rr.rates(1, irtiap)*rr.rates(index_rate, irw1);
       a(3) = -y(Cr48) * rr.rates(1, irw1) * rr.rates(index_rate, ircrgp);
       a(4) = -y(Cr48) * rr.rates(index_rate, irw1) * rr.rates(1, ircrgp);
       a(5) = -y(Cr48)*y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(1, irx1));
       a(6) =  y(Cr48)*y(He4) * rr.rates(1, ircrap)*rr.rates(index_rate, irx1);
       a(7) =  y(Fe52) * rr.rates(1, irx1) * rr.rates(index_rate, irfegp);
       a(8) =  y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(1, irfegp);

       dydt(Cr48) = dydt(Cr48) + esum8(a);

    }


    // fe52 reactions
    a(1) =  y(Cr48) * y(He4) * rr.rates(index_rate, ircrag);
    a(2) = -y(Fe52) * y(He4) * rr.rates(index_rate, irfeag);
    a(3) = -y(Fe52) * rr.rates(index_rate, irfega);
    a(4) =  y(Ni56) * rr.rates(index_rate, irniga);

    dydt(Fe52) = dydt(Fe52) + esum4(a);

    if (!deriva) {

       a(1) =  y(Cr48) * y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(index_rate, irx1));
       a(2) = -y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(index_rate, irfegp);
       a(3) = -y(Fe52) * y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(index_rate, iry1));
       a(4) =  y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(index_rate, irnigp);

       dydt(Fe52) = dydt(Fe52) + esum4(a);

    } else {

       a(1) =  y(Cr48)*y(He4) * rr.rates(index_rate, ircrap)*(1.0e0_rt-rr.rates(1, irx1));
       a(2) = -y(Cr48)*y(He4) * rr.rates(1, ircrap)*rr.rates(index_rate, irx1);
       a(3) = -y(Fe52) * rr.rates(1, irx1) * rr.rates(index_rate, irfegp);
       a(4) = -y(Fe52) * rr.rates(index_rate, irx1) * rr.rates(1, irfegp);
       a(5) = -y(Fe52)*y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(1, iry1));
       a(6) =  y(Fe52)*y(He4) * rr.rates(1, irfeap)*rr.rates(index_rate, iry1);
       a(7) =  y(Ni56) * rr.rates(1, iry1) * rr.rates(index_rate, irnigp);
       a(8) =  y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(1, irnigp);

       dydt(Fe52) = dydt(Fe52) + esum8(a);

    }


    // ni56 reactions
    a(1) =  y(Fe52) * y(He4) * rr.rates(index_rate, irfeag);
    a(2) = -y(Ni56) * rr.rates(index_rate, irniga);

    dydt(Ni56) = dydt(Ni56) + a(1) + a(2);

    if (!deriva) {

       a(1) =  y(Fe52) * y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(index_rate, iry1));
       a(2) = -y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(index_rate, irnigp);

       dydt(Ni56) = dydt(Ni56) + a(1) + a(2);

    } else {

       a(1) =  y(Fe52)*y(He4) * rr.rates(index_rate, irfeap)*(1.0e0_rt-rr.rates(1, iry1));
       a(2) = -y(Fe52)*y(He4) * rr.rates(1, irfeap)*rr.rates(index_rate, iry1);
       a(3) = -y(Ni56) * rr.rates(1, iry1) * rr.rates(index_rate, irnigp);
       a(4) = -y(Ni56) * rr.rates(index_rate, iry1) * rr.rates(1, irnigp);

       dydt(Ni56) = dydt(Ni56) + esum4(a);

    }

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

    RateEvaluator rr(rho, temp, y);

    // Call the RHS to actually get dydt.

    constexpr bool deriva = false;
    constexpr bool for_jacobian_tderiv = false;
    rhs(y, rr, ydot, deriva, for_jacobian_tderiv);

    // Instantaneous energy generation rate -- this needs molar fractions

    Real enuc;
    ener_gener_rate(ydot, enuc);

    // Get the neutrino losses

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Append the energy equation (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

    // Append the temperature equation

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

    RateEvaluator rr(rho, temp, y);

    // Species Jacobian elements with respect to other species

    dfdy_isotopes_aprox13(y, state, rr, jac);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc,j));
    }

    // Account for the thermal neutrino losses

    Real b1, sneut, dsneutdt, dsneutdd, snuda, snudz;

    sneut5(temp, rho, abar, zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    for (int j = 1; j <= NumSpec; ++j) {
       b1 = (-abar * abar * snuda + (zion[j-1] - zbar) * abar * snudz);
       jac.add(net_ienuc, j, -b1);
    }

    // Evaluate the Jacobian elements with respect to temperature by
    // calling the RHS using d(rate) / dT

    constexpr bool deriva = true;
    constexpr bool for_jacobian_tderiv = true;
    Array1D<Real, 1, neqs> yderivs;
    rhs(y, rr, yderivs, deriva, for_jacobian_tderiv);

    for (int i = 1; i <= NumSpec; ++i) {
        jac(i,net_itemp) = yderivs(i);
    }

    ener_gener_rate(yderivs, jac(net_ienuc,net_itemp));

    jac(net_ienuc,net_itemp) -= dsneutdt;

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
