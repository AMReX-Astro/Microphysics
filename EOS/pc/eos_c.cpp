#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

typedef double Real;

inline namespace literals {
    constexpr Real
    operator"" _rt( long double x )
    {
        return Real( x );
    }

    constexpr Real
    operator"" _rt( unsigned long long int x )
    {
        return Real( x );
    }
}

extern "C"
{
    // Inverse Fermi integral with q=1/2
    void ferinv7 (Real F, Real& X)
    {
        // Version 24.05.07
        // X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84, 101
        // Input: F
        // Output: X = X_q
        // Relative error:
        //        for X:    4.2e-9

        const Real A[6] = { 1.999266880833e4_rt,  5.702479099336e3_rt,  6.610132843877e2_rt,
                            3.818838129486e1_rt,  1.0_rt,               0.0_rt};
        const Real B[7] = { 1.771804140488e4_rt, -2.014785161019e3_rt,  9.130355392717e1_rt,
                           -1.670718177489e0_rt,  0.0_rt,               0.0_rt,
                            0.0_rt};
        const Real C[7] = {-1.277060388085e-2_rt,  7.187946804945e-2_rt, -4.262314235106e-1_rt,
                            4.997559426872e-1_rt, -1.285579118012e0_rt, -3.930805454272e-1_rt,
                            1.0_rt};
        const Real D[7] = {-9.745794806288e-3_rt,  5.485432756838e-2_rt, -3.29946624326e-1_rt,
                            4.077841975923e-1_rt, -1.145531476975e0_rt,  -6.067091689181e-2_rt,
                            0.0_rt};
        const int LA = 4;
        const int LB = 3;
        const int LD = 5;

        const int N = 1;

        if (F <= 0.0_rt) {
            printf("ferinv7: Non-positive argument\n");
            exit(1);
        }
        if (F < 4.0_rt) {
            Real T = F;
            Real UP = 0.0_rt;
            Real DOWN = 0.0_rt;
            for (int i = LA; i >= 0; --i) {
                UP = UP * T + A[i];
            }
            for (int i = LB; i >= 0; --i) {
                DOWN = DOWN * T + B[i];
            }
            X = std::log(T * UP / DOWN);
        }
        else {
            Real P = -1.0_rt / (0.5_rt + (Real) N); // = -1/(1+\nu) = power index
            Real T = std::pow(F, P); // t  - argument of the rational fraction
            Real UP = 0.0_rt;
            Real DOWN = 0.0_rt;
            for (int i = 6; i >= 0; --i) {
                UP = UP * T + C[i];
            }
            for (int i = LD; i >= 0; --i) {
                DOWN = DOWN * T + D[i];
            }
            Real R = UP / DOWN;
            X = R / T;
        }
    }

    void chemfit (Real DENS, Real TEMP, Real& CHI)
    {
        // Version 29.08.15
        // Fit to the chemical potential of free electron gas described in:
        //     G.Chabrier & A.Y.Potekhin, Phys.Rev.E 58, 4941 (1998)
        // Stems from CHEMFIT v.10.10.96. The main difference - derivatives.
        // Input:  DENS - electron density [a.u.=6.7483346e24 cm^{-3}],
        // TEMP - temperature [a.u.=2Ryd=3.1577e5 K]
        // Output: CHI = CMU1 / TEMR, where CMU1 = \mu-1 - chem.pot.w/o rest-energy

        const Real C13 = 1.0_rt / 3.0_rt;
        const Real PARA = 1.612_rt;
        const Real PARB = 6.192_rt;
        const Real PARC = 0.0944_rt;
        const Real PARF = 5.535_rt;
        const Real PARG = 0.698_rt;
        const Real XEPST = 228.0_rt; // the largest argument of e^{-X}

        Real DENR = DENS / 2.5733806e6_rt; // n_e in rel.un.=\lambda_{Compton}^{-3}
        Real TEMR = TEMP / 1.8778865e4_rt; // T in rel.un.=(mc^2/k)=5.93e9 K

        Real PF0 = std::pow(29.6088132_rt * DENR, C13); // Classical Fermi momentum
        Real TF;
        if (PF0 > 1.e-4_rt) {
            TF = std::sqrt(1.0_rt + PF0 * PF0) - 1.0_rt; // Fermi temperature
        }
        else {
            TF = 0.5_rt * PF0 * PF0;
        }

        Real THETA = TEMR / TF;
        Real THETA32 = THETA * std::sqrt(THETA);
        Real Q2 = 12.0_rt + 8.0_rt / THETA32;
        Real T1 = 0.0_rt;
        if (THETA < XEPST) {
            T1 = std::exp(-THETA);
        }
        Real U3 = T1 * T1 + PARA;
        Real THETAC = std::pow(THETA, PARC);
        Real THETAG = std::pow(THETA, PARG);
        Real D3 = PARB * THETAC * T1 * T1 + PARF * THETAG;
        Real Q3 = 1.365568127_rt - U3 / D3; // 1.365...=2/\pi^{1/3}
        Real Q1;
        if (THETA > 1.e-5_rt) {
            Q1 = 1.5_rt * T1 / (1.0_rt - T1);
        }
        else {
            Q1 = 1.5 / THETA;
        }
        Real SQT = std::sqrt(TEMR);
        Real G = (1.0_rt + Q2 * TEMR * Q3 + Q1 * SQT) * TEMR;
        Real H = (1.0_rt + 0.5 * TEMR / THETA) * (1.0_rt + Q2 * TEMR);
        Real CT = 1.0_rt + G / H;
        Real F = 2.0_rt * C13 / THETA32;
        Real X;
        ferinv7(F, X);
        CHI = X                      // Non-relativistic result
            - 1.5_rt * std::log(CT); // Relativistic fit
    }

    void blin9a (Real TEMP, Real CHI,
                 Real& W0, Real& W0DX, Real& W0DT, Real& W0DXX,
                 Real& W0DTT, Real& W0DXT,
                 Real& W1, Real& W1DX, Real& W1DT, Real& W1DXX,
                 Real& W1DTT, Real& W1DXT,
                 Real& W2, Real& W2DX, Real& W2DT, Real& W2DXX,
                 Real& W2DTT, Real& W2DXT,
                 Real& W0XXX, Real& W0XTT, Real& W0XXT)
    {
        // Version 19.01.10
        // First part of blin9: small CHI. Stems from blin9 v.24.12.08
        const Real AC[3][5] = {{0.37045057_rt,  0.41258437_rt,
                                9.777982e-2_rt, 5.3734153e-3_rt, 3.8746281e-5_rt},  // c_i^0
                               {0.39603109_rt,  0.69468795_rt,
                                0.22322760_rt,  1.5262934e-2_rt, 1.3081939e-4_rt},  // c_i^1
                               {0.76934619_rt,  1.7891437_rt,
                                0.70754974_rt,  5.6755672e-2_rt, 5.5571480e-4_rt}}; // c_i^2

        const Real AU[3][5] = {{0.43139881_rt,  1.7597537_rt,
                                4.10446540_rt,  7.7467038_rt, 13.457678_rt},  // \chi_i^0
                               {0.81763176_rt,  2.4723339_rt,
                                5.11600610_rt,  9.0441465_rt, 15.049882_rt},  // \chi_i^1
                               {1.25584610_rt,  3.2070406_rt,
                                6.12390820_rt, 10.3161260_rt, 16.597079_rt}}; // \chi_i^2

        const Real AA[3][5] = {{std::exp(-AU[0][0]), std::exp(-AU[0][1]),
                                std::exp(-AU[0][2]), std::exp(-AU[0][3]), std::exp(-AU[0][4])},  // \chi_i^0
                               {std::exp(-AU[1][0]), std::exp(-AU[1][1]),
                                std::exp(-AU[1][2]), std::exp(-AU[1][3]), std::exp(-AU[1][4])},  // \chi_i^1
                               {std::exp(-AU[2][0]), std::exp(-AU[2][1]),
                                std::exp(-AU[2][2]), std::exp(-AU[2][3]), std::exp(-AU[2][4])}}; // \chi_i^2

        for (int k = 0; k <= 2; ++k) {
            Real W = 0.0;
            Real WDX = 0.0;
            Real WDT = 0.0;
            Real WDXX = 0.0;
            Real WDTT = 0.0;
            Real WDXT = 0.0;
            Real WDXXX = 0.0;
            Real WDXTT = 0.0;
            Real WDXXT = 0.0;
            Real ECHI = std::exp(-CHI);

            for (int i = 0; i <= 4; ++i) {
                Real SQ = std::sqrt(1.0_rt + AU[k][i] * TEMP / 2.0_rt);
                Real DN = AA[k][i] + ECHI; // e^{-\chi_i}+e^{-\chi})

                W = W + AC[k][i] * SQ / DN;
                WDX = WDX + AC[k][i] * SQ / (DN * DN);
                WDT = WDT + AC[k][i] * AU[k][i] / (SQ * DN);
                WDXX = WDXX + AC[k][i] * SQ * (ECHI - AA[k][i]) / (DN * DN * DN);
                WDTT = WDTT - AC[k][i] * AU[k][i] * AU[k][i] / (DN * SQ * SQ * SQ);
                WDXT = WDXT + AC[k][i] * AU[k][i] / (SQ * DN * DN);
                WDXXX = WDXXX + AC[k][i] * SQ *
                        (ECHI * ECHI - 4.0_rt * ECHI * AA[k][i] + AA[k][i] * AA[k][i]) /
                        (DN * DN * DN * DN);
                WDXTT = WDXTT - AC[k][i] * AU[k][i] * AU[k][i] / (DN * DN * SQ * SQ * SQ);
                WDXXT = WDXXT + AC[k][i] * AU[k][i] * (ECHI - AA[k][i]) / (SQ * DN * DN * DN);
            }

            WDX = WDX * ECHI;
            WDT = WDT / 4.0_rt;
            WDXX = WDXX * ECHI;
            WDTT = WDTT / 16.0_rt;
            WDXT = WDXT / 4.0_rt * ECHI;
            WDXXX = WDXXX * ECHI;
            WDXTT = WDXTT * ECHI / 16.0_rt;
            WDXXT = WDXXT / 4.0_rt * ECHI;

            if (k == 0) {
                W0 = W;
                W0DX = WDX;
                W0DT = WDT;
                W0DXX = WDXX;
                W0DTT = WDTT;
                W0DXT = WDXT;
                W0XXX = WDXXX;
                W0XTT = WDXTT;
                W0XXT = WDXXT;
            }
            else if (k == 1) {
                W1 = W;
                W1DX = WDX;
                W1DT = WDT;
                W1DXX = WDXX;
                W1DTT = WDTT;
                W1DXT = WDXT;
            }
            else {
                W2 = W;
                W2DX = WDX;
                W2DT = WDT;
                W2DXX = WDXX;
                W2DTT = WDTT;
                W2DXT = WDXT;
            }
        }
    }

    void blin9b(Real TEMP, Real CHI,
                Real& W0, Real& W0DX, Real& W0DT, Real& W0DXX,
                Real& W0DTT, Real& W0DXT,
                Real& W1, Real& W1DX, Real& W1DT, Real& W1DXX,
                Real& W1DTT, Real& W1DXT,
                Real& W2, Real& W2DX, Real& W2DT, Real& W2DXX,
                Real& W2DTT, Real& W2DXT,
                Real& W0XXX, Real& W0XTT, Real& W0XXT)
    {
        // Version 19.01.10
        // Small syntax fix 15.03.13
        // Second part of BILN9: intermediate CHI. Stems from BLIN8 v.24.12.08
        const Real EPS = 1.e-3;

        const Real AX[5]  = {7.265351e-2_rt,  0.2694608_rt,
                             0.533122_rt,     0.7868801_rt,    0.9569313_rt}; // x_i
        const Real AXI[5] = {0.26356032_rt,   1.4134031_rt,
                             3.59642580_rt,   7.0858100_rt,   12.640801_rt}; // \xi_i
        const Real AH[5]  = {3.818735e-2_rt,  0.1256732_rt,
                             0.1986308_rt,    0.1976334_rt,    0.1065420_rt}; // H_i
        const Real AV[5]  = {0.29505869_rt,   0.32064856_rt,
                             7.3915570e-2_rt, 3.6087389e-3_rt, 2.3369894e-5_rt}; // \bar{V}_i

        if (CHI < EPS) {
            printf("BLIN9b: CHI is too small\n");
            exit(1);
        }

        for (int k = 0; k <= 2; ++k) {
            Real W = 0.0;
            Real WDX = 0.0;
            Real WDT = 0.0;
            Real WDXX = 0.0;
            Real WDTT = 0.0;
            Real WDXT = 0.0;
            Real WDXXX = 0.0;
            Real WDXTT = 0.0;
            Real WDXXT = 0.0;
            Real SQCHI = std::sqrt(CHI);

            for (int i = 0; i <= 4; ++i) {
                Real CE = AX[i] - 1.0_rt;
                Real ECHI = std::exp(CE * CHI);
                Real DE = 1.0_rt + ECHI;
                Real D = 1.0_rt + AX[i] * CHI * TEMP / 2.0_rt;
                Real H = std::pow(CHI, k + 1) * SQCHI * std::sqrt(D) / DE;
                Real HX = (k + 1.5_rt) / CHI + 0.25_rt * AX[i] * TEMP / D - ECHI * CE / DE;
                Real HDX = H * HX;
                Real HXX = (k + 1.5_rt) / (CHI * CHI) + 0.125_rt * (AX[i] * TEMP / D) * (AX[i] * TEMP / D) +
                           ECHI * (CE / DE) * (CE / DE);
                Real HDXX = HDX * HX - H * HXX;
                Real HT = 0.25_rt * AX[i] * CHI / D;
                Real HDT = H * HT;
                Real HDTT = -H * HT * HT;
                Real HTX = 1.0_rt / CHI - 0.5_rt * AX[i] * TEMP / D;
                Real HDXT = HDX * HT + HDT * HTX;
                Real HDXXT = HDXX * HT + HDX * HT * HTX + HDXT * HTX +
                             HDT * (0.25_rt * (AX[i] * TEMP / D) * (AX[i] * TEMP / D) -
                                    1.0_rt / (CHI * CHI));
                Real HDXTT = HDXT * HT - HDX * 0.125_rt * (AX[i] * CHI / D) * (AX[i] * CHI / D) + HDTT * HTX +
                             HDT * 0.5_rt * AX[i] * (TEMP * 0.5_rt * AX[i] * CHI / (D * D) - 1.0_rt / D);
                Real HXXX = (2 * k + 3) / (CHI * CHI * CHI) + 0.125_rt * (AX[i] * TEMP / D) * (AX[i] * TEMP / D) *
                            (AX[i] * TEMP / D) - ECHI * (1.0_rt - ECHI) * (CE / DE) * (CE / DE) * (CE / DE);
                Real HDXXX = HDXX * HX - 2.0_rt * HDX * HXX + H * HXXX;
                Real XICHI = AXI[i] + CHI;
                Real DXI = 1.0_rt + XICHI * TEMP / 2.0_rt;
                Real V = std::pow(XICHI, k) * std::sqrt(XICHI * DXI);
                Real VX= (k + 0.5_rt) / XICHI + 0.25_rt * TEMP / DXI;
                Real VDX = V * VX;
                Real VT = 0.25_rt * XICHI / DXI;
                Real VDT = V * VT;
                Real VXX = (k + 0.5_rt) / (XICHI * XICHI) + 0.125_rt * (TEMP / DXI) * (TEMP / DXI);
                Real VDXX = VDX * VX - V * VXX;
                Real VDXXX = VDXX * VX - 2.0_rt * VDX * VXX +
                             V * ((2 * k + 1) / (XICHI * XICHI * XICHI) +
                                  0.125_rt * (TEMP / DXI) * (TEMP / DXI) * (TEMP / DXI));
                Real VXXT = (1.0_rt - 0.5_rt * TEMP * XICHI / DXI) / DXI;
                Real VDTT = -V * VT * VT;
                Real VXT = 1.0_rt / XICHI - 0.5_rt * TEMP / DXI;
                Real VDXT = VDT * VXT + VDX * VT;
                Real VDXXT = VDXT * VX + VDX * 0.25_rt * VXXT - VDT * VXX - V * 0.25_rt * TEMP / DXI * VXXT;
                Real VDXTT = VDTT * VXT - VDT * 0.5_rt * VXXT + VDXT * VT -
                             VDX * 0.125_rt * (XICHI / DXI) * (XICHI / DXI);
                W = W + AH[i] * std::pow(AX[i], k) * H + AV[i] * V;
                WDX = WDX + AH[i] * std::pow(AX[i], k) * HDX + AV[i] * VDX;
                WDT = WDT + AH[i] * std::pow(AX[i], k) * HDT + AV[i] * VDT;
                WDXX = WDXX + AH[i] * std::pow(AX[i], k) * HDXX + AV[i] * VDXX;
                WDTT = WDTT + AH[i] * std::pow(AX[i], k) * HDTT + AV[i] * VDTT;
                WDXT = WDXT + AH[i] * std::pow(AX[i], k) * HDXT + AV[i] * VDXT;
                WDXXX = WDXXX + AH[i] * std::pow(AX[i], k) * HDXXX + AV[i] * VDXXX;
                WDXTT = WDXTT + AH[i] * std::pow(AX[i], k) * HDXTT + AV[i] * VDXTT;
                WDXXT = WDXXT + AH[i] * std::pow(AX[i], k) * HDXXT + AV[i] * VDXXT;
            }

            if (k == 0) {
                W0 = W;
                W0DX = WDX;
                W0DT = WDT;
                W0DXX = WDXX;
                W0DTT = WDTT;
                W0DXT = WDXT;
                W0XXX = WDXXX;
                W0XTT = WDXTT;
                W0XXT = WDXXT;
            }
            else if (k == 1) {
                W1 = W;
                W1DX = WDX;
                W1DT = WDT;
                W1DXX = WDXX;
                W1DTT = WDTT;
                W1DXT = WDXT;
            }
            else {
                W2 = W;
                W2DX = WDX;
                W2DT = WDT;
                W2DXX = WDXX;
                W2DTT = WDTT;
                W2DXT = WDXT;
            }
        }
    }

    void blin9c (Real TEMP, Real CHI,
                 Real& W0, Real& W0DX, Real& W0DT, Real& W0DXX, Real& W0DTT, Real& W0DXT,
                 Real& W1, Real& W1DX, Real& W1DT, Real& W1DXX, Real& W1DTT, Real& W1DXT,
                 Real& W2, Real& W2DX, Real& W2DT, Real& W2DXX, Real& W2DTT, Real& W2DXT,
                 Real& W0XXX, Real& W0XTT, Real& W0XXT)
    {
        // Version 19.01.10
        // Third part of BILN9: large CHI. Stems from BLIN8 v.24.12.08
        const Real PI = 3.141592653_rt;
        const Real PI26 = PI * PI / 6.0;

        Real AM[3], AMDX[3], AMDT[3], AMDXX[3], AMDTT[3], AMDXT[3];

        if (CHI * TEMP < 0.1_rt) {

            for (int k = 0; k <= 2; ++k) {
                Real W = 0.0_rt;
                Real WDX = 0.0_rt;
                Real WDT = 0.0_rt;
                Real WDXX = 0.0_rt;
                Real WDTT = 0.0_rt;
                Real WDXT = 0.0_rt;
                Real WDXXX = 0.0_rt;
                Real WDXTT = 0.0_rt;
                Real WDXXT = 0.0_rt;

                Real C;

                for (int j = 0; j <= 4; ++j) { // for nonrel.Fermi integrals from k+1/2 to k+4.5
                    Real CNU = k + j + 0.5_rt; // nonrelativistic Fermi integral index \nu
                    Real CHINU = std::pow(CHI, k + j) * std::sqrt(CHI); // \chi^\nu
                    Real F = CHINU * (CHI / (CNU + 1.0_rt) + PI26 * CNU / CHI + // nonrel.Fermi
                                      0.7_rt * PI26 * PI26 * CNU * (CNU - 1.0_rt) *
                                      (CNU - 2.0_rt) / (CHI * CHI * CHI));
                    Real FDX = CHINU * (1.0_rt + PI26 * CNU * (CNU - 1.0_rt) / (CHI * CHI) +
                                        0.7_rt * PI26 * PI26 * CNU * (CNU - 1.0_rt) * (CNU - 2.0_rt)
                                        * (CNU - 3.0_rt) / (CHI * CHI * CHI * CHI));
                    Real FDXX = CHINU / CHI * CNU *
                                (1.0_rt + PI26 * (CNU - 1.0_rt) *
                                 (CNU - 2.0_rt) / (CHI * CHI) +
                                 0.7_rt * PI26 * PI26 * (CNU - 1.0_rt) * (CNU - 2.0_rt) *
                                 (CNU - 3.0_rt) * (CNU - 4.0_rt) / (CHI * CHI * CHI * CHI));
                    Real FDXXX = CHINU / (CHI * CHI) * CNU * (CNU - 1.0_rt) *
                                 (1.0_rt + PI26 * (CNU - 2.0_rt) * (CNU - 3.0_rt) / (CHI * CHI) +
                                  0.7_rt * PI26 * PI26 * (CNU - 2.0_rt) * (CNU - 3.0_rt) *
                                  (CNU - 4.0_rt) * (CNU - 5.0_rt) / (CHI * CHI * CHI * CHI));

                    if (j == 0) {
                        W = F;
                        WDX = FDX;
                        WDXX = FDXX;
                        WDXXX = FDXXX;
                    }
                    else if (j == 1) {
                        C = 0.25_rt * TEMP;
                        W = W + C * F; // Fermi-Dirac, expressed through Fermi
                        WDX = WDX + C * FDX;
                        WDXX = WDXX + C * FDXX;
                        WDT = F / 4.0_rt;
                        WDXT = FDX / 4.0_rt;
                        WDTT = 0.0_rt;
                        WDXXX = WDXXX + C * FDXXX;
                        WDXXT = FDXX / 4.0_rt;
                        WDXTT = 0.0_rt;
                    }
                    else {
                        C = -C / j * (2 * j - 3) / 4.0_rt * TEMP;
                        W = W + C * F;
                        WDX = WDX + C * FDX;
                        WDT = WDT + C * j / TEMP * F;
                        WDXX = WDXX + C * FDXX;
                        WDTT = WDTT + C * j * (j - 1) / (TEMP * TEMP) * F;
                        WDXT = WDXT + C * j / TEMP * FDX;
                        WDXXX = WDXXX + C * FDXXX;
                        WDXTT = WDXTT + C * j * (j - 1) / (TEMP * TEMP) * FDX;
                        WDXXT = WDXXT + C * j / TEMP * FDXX;
                    }
                }

                if (k == 0) {
                    W0 = W;
                    W0DX = WDX;
                    W0DT = WDT;
                    W0DXX = WDXX;
                    W0DTT = WDTT;
                    W0DXT = WDXT;
                    W0XXX = WDXXX;
                    W0XTT = WDXTT;
                    W0XXT = WDXXT;
                }
                else if (k == 1) {
                    W1 = W;
                    W1DX = WDX;
                    W1DT = WDT;
                    W1DXX = WDXX;
                    W1DTT = WDTT;
                    W1DXT = WDXT;
                }
                else {
                    W2 = W;
                    W2DX = WDX;
                    W2DT = WDT;
                    W2DXX = WDXX;
                    W2DTT = WDTT;
                    W2DXT = WDXT;
                }
            }

        }
        else { // CHI > 14, CHI * TEMP > 0.1: general high-\chi expansion

            Real D = 1.0_rt + CHI * TEMP / 2.0_rt;
            Real R = std::sqrt(CHI * D);
            Real RX = 0.5_rt / CHI + 0.25_rt * TEMP / D;
            Real RDX = R * RX;
            Real RDT = 0.25_rt * CHI * CHI / R;
            Real RXX = -0.5_rt / (CHI * CHI) - 0.125_rt * (TEMP / D) * (TEMP / D);
            Real RDXX = RDX * RX + R * RXX;
            Real RDTT = -0.25_rt * RDT * CHI / D;
            Real RXT = 0.25_rt / D - 0.125_rt * CHI * TEMP / (D * D);
            Real RDXT = RDT * RX + R * RXT;
            Real RXXX = 1.0_rt / (CHI * CHI * CHI) + 0.125_rt * (TEMP / D) * (TEMP / D) * (TEMP / D);
            Real RDXXX = RDXX * RX + 2.0_rt * RDX * RXX + R * RXXX;
            Real RXTT = -0.25_rt / (D * D) * CHI + 0.125_rt * CHI * CHI * TEMP / (D * D * D);
            Real RDXTT = RDTT * RX + 2.0_rt * RDT * RXT + R * RXTT;
            Real RXXT = -RXT * TEMP / D;
            Real RDXXT = RDXT * RX + RDX * RXT + RDT * RXX + R * RXXT;

            Real AMDXXX, AMDXTT, AMDXXT;

            for (int k = 0; k <= 2; ++k) {
                Real DM = k + 0.5_rt + (k + 1.0_rt) * CHI * TEMP / 2.0_rt;
                AM[k] = std::pow(CHI, k) * DM / R;
                Real FMX1 = 0.5_rt * (k + 1.0_rt) * TEMP / DM;
                Real FMX2 = 0.25_rt * TEMP / D;
                Real FMX = (k - 0.5_rt) / CHI + FMX1 - FMX2;
                AMDX[k] = AM[k] * FMX;
                Real CkM = 0.5_rt * (k + 1.0_rt) / DM;
                Real FMT1 = CkM * CHI;
                Real FMT2 = 0.25_rt * CHI / D;
                Real FMT = FMT1 - FMT2;
                AMDT[k] = AM[k] * FMT;
                Real FMXX = -(k - 0.5_rt) / (CHI * CHI) - FMX1 * FMX1 + 2.0_rt * FMX2 * FMX2;
                AMDXX[k] = AMDX[k] * FMX + AM[k] * FMXX;
                Real FMTT = 2.0_rt * FMT2 * FMT2 - FMT1 * FMT1;
                AMDTT[k] = AMDT[k] * FMT + AM[k] * FMTT;
                AMDXT[k] = AMDX[k] * FMT + AM[k] * (CkM * (1.0_rt - CkM * CHI * TEMP) -
                                                    0.25_rt / D + 0.125_rt * CHI * TEMP / (D * D));

                if (k == 0) {
                    Real FMXXX = (2 * k - 1) / (CHI * CHI * CHI) + 2.0_rt * FMX1 * FMX1 * FMX1 -
                                 8.0_rt * FMX2 * FMX2 * FMX2;
                    AMDXXX = AMDXX[k] * FMX + 2.0_rt * AMDX[k] * FMXX + AM[k] * FMXXX;
                    Real FMT1DX = CkM - TEMP * CHI * CkM * CkM;
                    Real FMT2DX = (0.25_rt - CHI * TEMP * 0.125_rt / D) / D;
                    Real FMXT = FMT1DX - FMT2DX;
                    Real FMTTX = 4.0_rt * FMT2 * FMT2DX - 2.0_rt * FMT1 * FMT1DX;
                    AMDXTT = AMDXT[k] * FMT + AMDT[k] * FMXT + AMDX[k] * FMTT + AM[k] * FMTTX;
                    Real FMX1DT = CkM - CHI * TEMP * CkM * CkM;
                    Real FMX2DT = 0.25_rt / D * (1.0_rt - 0.5_rt * CHI * TEMP / D);
                    Real FMXXT = 4.0_rt * FMX2 * FMX2DT - 2.0_rt * FMX1 * FMX1DT;
                    AMDXXT = AMDXT[k] * FMX + AMDX[k] * FMXT + AMDT[k] * FMXX + AM[k] * FMXXT;
                }
            }

            Real SQ2T = std::sqrt(2.0_rt * TEMP);
            Real A = 1.0_rt + CHI * TEMP + SQ2T * R;
            Real ADX = TEMP + SQ2T * RDX;
            Real ADT = CHI + R / SQ2T + SQ2T * RDT;
            Real ADXX = SQ2T * RDXX;
            Real ADTT = -R / (SQ2T * SQ2T * SQ2T) + 2.0_rt / SQ2T * RDT + SQ2T * RDTT;
            Real ADXT = 1.0_rt + RDX / SQ2T + SQ2T * RDXT;
            Real ADXTT = -RDX / (SQ2T * SQ2T * SQ2T) + 2.0_rt / SQ2T * RDXT + SQ2T * RDXTT;
            Real ADXXT = RDXX / SQ2T + SQ2T * RDXXT;
            Real XT1 = CHI + 1.0_rt / TEMP;
            Real Aln = std::log(A);
            Real FJ0 = 0.5_rt * XT1 * R - Aln / (SQ2T * SQ2T * SQ2T);
            Real ASQ3 = A * SQ2T * SQ2T * SQ2T;
            Real ASQ3DX = ADX * SQ2T * SQ2T * SQ2T;
            Real FJ0DX = 0.5_rt * (R + XT1 * RDX) - ADX / ASQ3;
            Real FJ0DT = 0.5_rt * (XT1 * RDT - R / (TEMP * TEMP)) - ADT / ASQ3 +
                         0.75_rt / (TEMP * TEMP * SQ2T) * Aln;
            Real FJ0DXX = RDX + 0.5_rt * XT1 * RDXX + (ADX / A) * (ADX / A) / (SQ2T * SQ2T * SQ2T) - ADXX / ASQ3;
            Real FJ0DTT = R / (TEMP * TEMP * TEMP) - RDT / (TEMP * TEMP) + 0.5_rt * XT1 * RDTT +
                          3.0_rt / (ASQ3 * TEMP) * ADT +
                          (ADT / A) * (ADT / A) / (SQ2T * SQ2T * SQ2T) - ADTT / ASQ3 -
                          1.875_rt / (TEMP * TEMP * TEMP * SQ2T) * Aln;
            Real BXT = 1.5_rt / TEMP * ADX + ADX * ADT / A - ADXT;
            Real BXXT = 1.5_rt / TEMP * ADXX + (ADXX * ADT + ADX * ADXT) / A -
                        (ADX / A) * (ADX / A) * ADT - ADXXT;
            Real FJ0DXT = 0.5_rt * (RDT - RDX / (TEMP * TEMP) + XT1 * RDXT) + BXT / ASQ3;
            Real FJ0XXX = RDXX * 1.5_rt + 0.5_rt * XT1 * RDXXX +
                          (2.0_rt * ADX * (ADXX / A - (ADX / A) * (ADX / A)) -
                          SQ2T * RDXXX + ADXX / ASQ3 * ASQ3DX) / ASQ3;
            Real FJ0XTT = RDX / (TEMP * TEMP * TEMP) - RDXT / (TEMP * TEMP) + 0.5_rt * (RDTT + XT1 * RDXTT) +
                          3.0_rt / TEMP * (ADXT - ADT / ASQ3 * ASQ3DX) / ASQ3 +
                          (2.0_rt * ADT * (ADXT / A - ADT * ADX / (A * A)) -
                           ADXTT + ADTT * ASQ3DX / ASQ3) / ASQ3 - 1.875_rt / (TEMP * TEMP * TEMP * SQ2T) * ADX / A;
            Real FJ0XXT = 0.5_rt * (RDXT - RDXX / (TEMP * TEMP) + RDXT + XT1 * RDXXT) +
                          (BXXT - BXT * ASQ3DX / ASQ3) / ASQ3;

            W0 = FJ0 + PI26 * AM[0];
            W0DX = FJ0DX + PI26 * AMDX[0];
            W0DT = FJ0DT + PI26 * AMDT[0];
            W0DXX = FJ0DXX + PI26 * AMDXX[0];
            W0DTT = FJ0DTT + PI26 * AMDTT[0];
            W0DXT = FJ0DXT + PI26 * AMDXT[0];
            W0XXX = FJ0XXX + PI26 * AMDXXX;
            W0XTT = FJ0XTT + PI26 * AMDXTT;
            W0XXT = FJ0XXT + PI26 * AMDXXT;

            Real FJ1 = (R * R * R / 1.5_rt - FJ0) / TEMP;
            Real FJ1DX = (2.0_rt * R * R * RDX - FJ0DX) / TEMP;
            Real FJ1DT = (2.0_rt * R * R * RDT - FJ0DT - FJ1) / TEMP;
            Real FJ1DXX = (4.0_rt * R * RDX * RDX + 2.0_rt * R * R * RDXX - FJ0DXX) / TEMP;
            Real FJ1DTT = (4.0_rt * R * RDT * RDX + 2.0_rt * R * R * RDTT - FJ0DTT - 2.0_rt * FJ1DT) / TEMP;
            Real FJ1DXT = (4.0_rt * R * RDX * RDT + 2.0_rt * R * R * RDXT - FJ0DXT - FJ1DX) / TEMP;

            W1 = FJ1 + PI26 * AM[1];
            W1DX = FJ1DX + PI26 * AMDX[1];
            W1DT = FJ1DT + PI26 * AMDT[1];
            W1DXX = FJ1DXX + PI26 * AMDXX[1];
            W1DTT = FJ1DTT + PI26 * AMDTT[1];
            W1DXT = FJ1DXT + PI26 * AMDXT[1];

            Real FJ2 = (0.5_rt * CHI * R * R * R - 1.25_rt * FJ1) / TEMP;
            Real FJ2DX = (0.5_rt * R * R * R + 1.5_rt * CHI * R * R * RDX - 1.25_rt * FJ1DX) / TEMP;
            Real FJ2DT = (1.5_rt * CHI * R * R * RDT - 1.25_rt * FJ1DT - FJ2) / TEMP;
            Real FJ2DXX = (3.0_rt * R * RDX * (R + CHI * RDX) + 1.5_rt * CHI * R * R * RDXX -
                          1.25_rt * FJ1DXX) / TEMP;
            Real FJ2DTT = (3.0_rt * CHI * R * (RDT * RDT + 0.5_rt * R * RDTT) -
                          1.25_rt * FJ1DTT - 2.0_rt * FJ2DT) / TEMP;
            Real FJ2DXT = (1.5_rt * R * RDT * (R + 2.0_rt * CHI * RDX) + 1.5_rt * CHI * R * R * RDXT -
                           1.25_rt * FJ1DXT - FJ2DX) / TEMP;

            W2 = FJ2 + PI26 * AM[2];
            W2DX = FJ2DX + PI26 * AMDX[2];
            W2DT = FJ2DT + PI26 * AMDT[2];
            W2DXX = FJ2DXX + PI26 * AMDXX[2];
            W2DTT = FJ2DTT + PI26 * AMDTT[2];
            W2DXT = FJ2DXT + PI26 * AMDXT[2];
        }

    }

    void fermi10 (Real X, Real XMAX, Real& FP, Real& FM)
    {
        // Version 20.01.10
        // Fermi distribution function and its 3 derivatives
        // Input: X - argument f(x)
        //        XMAX - max|X| where it is assumed that 0 < f(x) < 1.
        // Output: FP = f(x)
        //         FM = 1-f(x)
        if (X > XMAX) {
            FP = 0.0_rt;
            FM = 1.0_rt;
        }
        else if (X < -XMAX) {
            FP = 1.0_rt;
            FM = 0.0_rt;
        }
        else {
            FP = 1.0 / (std::exp(X) + 1.0_rt);
            FM = 1.0 - FP;
        }
    }

    void blin9 (Real TEMP, Real CHI,
                Real& W0, Real& W0DX, Real& W0DT, Real& W0DXX, Real& W0DTT, Real& W0DXT,
                Real& W1, Real& W1DX, Real& W1DT, Real& W1DXX, Real& W1DTT, Real& W1DXT,
                Real& W2, Real& W2DX, Real& W2DT, Real& W2DXX, Real& W2DTT, Real& W2DXT,
                Real& W0XXX, Real& W0XTT, Real& W0XXT)
    {
        // Version 21.01.10
        // Stems from BLIN8 v.24.12.08
        // Difference - smooth matching of different CHI ranges
        // Input: TEMP=T/mc^2; CHI=(\mu-mc^2)/T
        // Output: Wk - Fermi-Dirac integral of the order k+1/2
        //         WkDX=dWk/dCHI, WkDT = dWk/dT, WkDXX=d^2 Wk / d CHI^2,
        //         WkDTT=d^2 Wk / d T^2, WkDXT=d^2 Wk /dCHIdT,
        //         W0XXX=d^3 W0 / d CHI^3, W0XTT=d^3 W0 /(d CHI d^2 T),
        //         W0XXT=d^3 W0 /dCHI^2 dT

        const Real CHI1 = 0.6_rt;
        const Real CHI2 = 14.0_rt;
        const Real XMAX = 30.0_rt;
        const Real DCHI1 = 0.1_rt;
        const Real DCHI2 = CHI2 - CHI1 - DCHI1;
        const Real XSCAL1 = XMAX / DCHI1;
        const Real XSCAL2 = XMAX / DCHI2;

        Real X1 = (CHI - CHI1) * XSCAL1;
        Real X2 = (CHI - CHI2) * XSCAL2;

        if (X1 < - XMAX) {

            blin9a(TEMP, CHI,
                   W0, W0DX, W0DT, W0DXX, W0DTT, W0DXT,
                   W1, W1DX, W1DT, W1DXX, W1DTT, W1DXT,
                   W2, W2DX, W2DT, W2DXX, W2DTT, W2DXT,
                   W0XXX, W0XTT, W0XXT);

        }
        else if (X2 < XMAX) { // match two fits

            Real W0a, W0DXa, W0DTa, W0DXXa, W0DTTa, W0DXTa,
                 W1a, W1DXa, W1DTa, W1DXXa, W1DTTa, W1DXTa,
                 W2a, W2DXa, W2DTa, W2DXXa, W2DTTa, W2DXTa,
                 W0XXXa, W0XTTa, W0XXTa;

            Real W0b, W0DXb, W0DTb, W0DXXb, W0DTTb, W0DXTb,
                 W1b, W1DXb, W1DTb, W1DXXb, W1DTTb, W1DXTb,
                 W2b, W2DXb, W2DTb, W2DXXb, W2DTTb, W2DXTb,
                 W0XXXb, W0XTTb, W0XXTb;

            Real FP, FM;

            if (X1 < XMAX) { // match fits "a" and "b"

                fermi10(X1, XMAX, FP, FM);
                blin9a(TEMP, CHI,
                       W0a, W0DXa, W0DTa, W0DXXa, W0DTTa, W0DXTa,
                       W1a, W1DXa, W1DTa, W1DXXa, W1DTTa, W1DXTa,
                       W2a, W2DXa, W2DTa, W2DXXa, W2DTTa, W2DXTa,
                       W0XXXa, W0XTTa, W0XXTa);
                blin9b(TEMP, CHI,
                       W0b, W0DXb, W0DTb, W0DXXb, W0DTTb, W0DXTb,
                       W1b, W1DXb, W1DTb, W1DXXb, W1DTTb, W1DXTb,
                       W2b, W2DXb, W2DTb, W2DXXb, W2DTTb, W2DXTb,
                       W0XXXb, W0XTTb, W0XXTb);

            }
            else { // match fits "b" and "c"

                fermi10(X2, XMAX, FP, FM);
                blin9b(TEMP, CHI,
                       W0a, W0DXa, W0DTa, W0DXXa, W0DTTa, W0DXTa,
                       W1a, W1DXa, W1DTa, W1DXXa, W1DTTa, W1DXTa,
                       W2a, W2DXa, W2DTa, W2DXXa, W2DTTa, W2DXTa,
                       W0XXXa, W0XTTa, W0XXTa);
                blin9c(TEMP, CHI,
                       W0b, W0DXb, W0DTb, W0DXXb, W0DTTb, W0DXTb,
                       W1b, W1DXb, W1DTb, W1DXXb, W1DTTb, W1DXTb,
                       W2b, W2DXb, W2DTb, W2DXXb, W2DTTb, W2DXTb,
                       W0XXXb, W0XTTb, W0XXTb);

            }

            W0 = W0a * FP + W0b * FM;
            W0DX = W0DXa * FP + W0DXb * FM;
            W0DT = W0DTa * FP + W0DTb * FM;
            W0DXX = W0DXXa * FP + W0DXXb * FM;
            W0DTT = W0DTTa * FP + W0DTTb * FM;
            W0DXT = W0DXTa * FP + W0DXTb * FM;
            W0XXX = W0XXXa * FP + W0XXXb * FM;
            W0XTT = W0XTTa * FP + W0XTTb * FM;
            W0XXT = W0XXTa * FP + W0XXTb * FM;
            W1 = W1a * FP + W1b * FM;
            W1DX = W1DXa * FP + W1DXb * FM;
            W1DT = W1DTa * FP + W1DTb * FM;
            W1DXX = W1DXXa * FP + W1DXXb * FM;
            W1DTT = W1DTTa * FP + W1DTTb * FM;
            W1DXT = W1DXTa * FP + W1DXTb * FM;
            W2 = W2a * FP + W2b * FM;
            W2DX = W2DXa * FP + W2DXb * FM;
            W2DT = W2DTa * FP + W2DTb * FM;
            W2DXX = W2DXXa * FP + W2DXXb * FM;
            W2DTT = W2DTTa * FP + W2DTTb * FM;
            W2DXT = W2DXTa * FP + W2DXTb * FM;

        }
        else {

            blin9c(TEMP, CHI,
                   W0, W0DX, W0DT, W0DXX, W0DTT, W0DXT,
                   W1, W1DX, W1DT, W1DXX, W1DTT, W1DXT,
                   W2, W2DX, W2DT, W2DXX, W2DTT, W2DXT,
                   W0XXX, W0XTT, W0XXT);

        }
    }
}
