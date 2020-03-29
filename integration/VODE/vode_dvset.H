#ifndef _cuvode_dvset_H_
#define _cuvode_dvset_H_

#include <vode_type.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void dvset (dvode_t& vstate)
{

    // dvset is called by dvstep and sets coefficients for use there.
    //
    // For each order NQ, the coefficients in EL are calculated by use of
    // the generating polynomial lambda(x), with coefficients EL(i).
    //      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
    // Since we are using the BDF method, this is
    //                                     NQ-1
    //     lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
    //                                     i = 1
    // xi(i) is defined by
    // H*xi(i) = t sub n  -  t sub (n-i)
    //         = H + TAU(1) + TAU(2) + ... TAU(i-1).

    const Real CORTES = 0.1e0_rt;

    Real FLOTL = vstate.L;
    int NQM1 = vstate.NQ - 1;
    int NQM2 = vstate.NQ - 2;

    for (int i = 3; i <= vstate.L; ++i) {
        vstate.el(i) = 0.0_rt;
    }

    vstate.el(1) = 1.0_rt;
    vstate.el(2) = 1.0_rt;
    Real ALPH0 = -1.0_rt;
    Real AHATN0 = -1.0_rt;
    Real HSUM = vstate.H;
    Real RXI = 1.0_rt;
    Real RXIS = 1.0_rt;

    if (vstate.NQ != 1) {
        for (int j = 1; j <= NQM2; ++j) {
            // In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)).
            HSUM += vstate.tau(j);
            RXI = vstate.H / HSUM;
            ALPH0 -= 1.0_rt / (j+1);
            for (int iback = 1; iback <= j+1; ++iback) {
                int i = (j + 3) - iback;
                vstate.el(i) += vstate.el(i-1) * RXI;
            }
        }

        ALPH0 -= 1.0_rt / vstate.NQ;
        RXIS = -vstate.el(2) - ALPH0;
        HSUM += vstate.tau(NQM1);
        RXI = vstate.H / HSUM;
        AHATN0 = -vstate.el(2) - RXI;
        for (int iback = 1; iback <= vstate.NQ; ++iback) {
            int i = (vstate.NQ + 2) - iback;
            vstate.el(i) += vstate.el(i-1) * RXIS;
        }
    }

    Real T1 = 1.0_rt - AHATN0 + ALPH0;
    Real T2 = 1.0_rt + vstate.NQ * T1;
    vstate.tq(2) = std::abs(ALPH0 * T2 / T1);
    vstate.tq(5) = std::abs(T2 / (vstate.el(vstate.L) * RXI / RXIS));

    if (vstate.NQWAIT == 1) {
        Real CNQM1 = RXIS / vstate.el(vstate.L);
        Real T3 = ALPH0 + 1.0_rt / vstate.NQ;
        Real T4 = AHATN0 + RXI;
        Real ELP = T3 / (1.0_rt - T4 + T3);
        vstate.tq(1) = std::abs(ELP / CNQM1);
        HSUM += vstate.tau(vstate.NQ);
        RXI = vstate.H / HSUM;
        Real T5 = ALPH0 - 1.0_rt / (vstate.NQ+1);
        Real T6 = AHATN0 - RXI;
        ELP = T2 / (1.0_rt - T6 + T5);
        vstate.tq(3) = std::abs(ELP * RXI * (FLOTL + 1.0_rt) * T5);
    }

    vstate.tq(4) = CORTES * vstate.tq(2);

}

#endif
