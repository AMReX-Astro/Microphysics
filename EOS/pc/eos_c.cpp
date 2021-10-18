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

    void excor7 (double RS, double GAME,
                 double& FXC, double& UXC, double& PXC,
                 double& CVXC, double& SXC, double& PDTXC,
                 double& PDRXC)
    {
        // Version 09.06.07
        // Accuracy-loss cut-off modified on 10.08.16
        // Exchange-correlation contribution for the electron gas
        // Stems from TANAKA1 v.03.03.96. Added derivatives.
        // Input: RS - electron density parameter =electron-sphere radius in a.u.
        //        GAME - electron Coulomb coupling parameter
        // Output: FXC - excess free energy of e-liquid per kT per one electron
        //               according to Tanaka & Ichimaru 85-87 and Ichimaru 93
        //         UXC - internal energy contr.[per 1 electron, kT]
        //         PXC - pressure contribution divided by (n_e kT)
        //         CVXC - heat capacity divided by N_e k
        //         SXC - entropy divided by N_e k
        //         PDTXC,PDRXC = PXC + d PXC / d ln(T,\rho)
        const Real EPS = 1.e-8_rt; // 10.08.16

        Real THETA = 0.543_rt * RS / GAME; // non-relativistic degeneracy parameter
        Real SQTH = std::sqrt(THETA);
        Real THETA2 = THETA * THETA;
        Real THETA3 = THETA2 * THETA;
        Real THETA4 = THETA3 * THETA;

        Real T1, T1DH, T1DHH, T2, T2DH, T2DHH;

        if (THETA > .005_rt) {
            Real CHT1 = std::cosh(1.0_rt / THETA);
            Real SHT1 = std::sinh(1.0_rt / THETA);
            Real CHT2 = std::cosh(1.0_rt / SQTH);
            Real SHT2 = std::sinh(1.0_rt / SQTH);
            T1 = SHT1 / CHT1; // tanh(1.0_rt / THETA)
            T2 = SHT2 / CHT2; // tanh(1.0_rt / sqrt(THETA))
            T1DH = -1.0_rt / ((THETA * CHT1) * (THETA * CHT1)); // d T1 / d\theta
            T1DHH = 2.0_rt / ((THETA * CHT1) * (THETA * CHT1) * (THETA * CHT1)) *
                    (CHT1 - SHT1 / THETA);
            T2DH =  -0.5_rt * SQTH / ((THETA * CHT2) * (THETA * CHT2));
            T2DHH = (0.75_rt * SQTH * CHT2 - 0.5_rt * SHT2) /
                    ((THETA * CHT2) * (THETA * CHT2) * (THETA * CHT2));
        }
        else {
            T1 = 1.0_rt;
            T2 = 1.0_rt;
            T1DH = 0.0_rt;
            T2DH = 0.0_rt;
            T1DHH = 0.0_rt;
            T2DHH = 0.0_rt;
        }

        Real A0 = 0.75_rt + 3.04363_rt * THETA2 - 0.09227_rt * THETA3 + 1.7035_rt * THETA4;
        Real A0DH = 6.08726_rt * THETA - 0.27681_rt * THETA2 + 6.814_rt * THETA3;
        Real A0DHH = 6.08726_rt - 0.55362_rt * THETA + 20.442_rt * THETA2;
        Real A1 = 1.0_rt + 8.31051_rt * THETA2 + 5.1105_rt * THETA4;
        Real A1DH = 16.62102_rt * THETA + 20.442_rt * THETA3;
        Real A1DHH = 16.62102_rt + 61.326_rt * THETA2;
        Real A = 0.610887_rt * A0 / A1 * T1; // HF fit of Perrot and Dharma - wardana
        Real AH = A0DH / A0 - A1DH / A1 + T1DH / T1;
        Real ADH = A * AH;
        Real ADHH = ADH * AH + A * (A0DHH / A0 - (A0DH / A0) * (A0DH / A0) -
                                    A1DHH / A1 + (A1DH / A1) * (A1DH / A1) +
                                    T1DHH / T1 - (T1DH / T1) * (T1DH / T1));
        Real B0 = 0.341308_rt + 12.070873_rt * THETA2 + 1.148889_rt * THETA4;
        Real B0DH = 24.141746_rt * THETA + 4.595556_rt * THETA3;
        Real B0DHH = 24.141746_rt + 13.786668_rt * THETA2;
        Real B1 = 1.0_rt + 10.495346_rt * THETA2 + 1.326623 * THETA4;
        Real B1DH = 20.990692_rt * THETA + 5.306492 * THETA3;
        Real B1DHH = 20.990692_rt + 15.919476_rt * THETA2;
        Real B = SQTH * T2 * B0 / B1;
        Real BH = 0.5_rt / THETA + T2DH / T2 + B0DH / B0 - B1DH / B1;
        Real BDH = B * BH;
        Real BDHH = BDH * BH + B * (-0.5_rt / THETA2 + T2DHH / T2 - (T2DH / T2) * (T2DH / T2) +
                                    B0DHH / B0 - (B0DH / B0) * (B0DH / B0) - B1DHH / B1 +
                                    (B1DH / B1) * (B1DH / B1));
        Real D0 = 0.614925_rt + 16.996055_rt * THETA2 + 1.489056_rt * THETA4;
        Real D0DH = 33.99211_rt * THETA + 5.956224_rt * THETA3;
        Real D0DHH = 33.99211_rt + 17.868672_rt * THETA2;
        Real D1 = 1.0_rt + 10.10935_rt * THETA2 + 1.22184_rt * THETA4;
        Real D1DH = 20.2187_rt * THETA + 4.88736_rt * THETA3;
        Real D1DHH = 20.2187_rt + 14.66208_rt * THETA2;
        Real D = SQTH * T2 * D0 / D1;
        Real DH = 0.5_rt / THETA + T2DH / T2 + D0DH / D0 - D1DH / D1;
        Real DDH = D * DH;
        Real DDHH = DDH * DH + D * (-0.5_rt / THETA2 + T2DHH / T2 - (T2DH / T2) * (T2DH / T2) +
                                    D0DHH / D0 - (D0DH / D0) * (D0DH / D0) - D1DHH / D1 +
                                    (D1DH / D1) * (D1DH / D1));
        Real E0 = 0.539409_rt + 2.522206_rt * THETA2 + 0.178484_rt * THETA4;
        Real E0DH = 5.044412_rt * THETA + 0.713936_rt * THETA3;
        Real E0DHH = 5.044412_rt + 2.141808_rt * THETA2;
        Real E1 = 1.0_rt + 2.555501_rt * THETA2 + 0.146319_rt * THETA4;
        Real E1DH = 5.111002_rt * THETA + 0.585276_rt * THETA3;
        Real E1DHH = 5.111002_rt + 1.755828_rt * THETA2;
        Real E = THETA * T1 * E0 / E1;
        Real EH = 1.0_rt / THETA + T1DH / T1 + E0DH / E0 - E1DH / E1;
        Real EDH = E * EH;
        Real EDHH = EDH * EH + E * (T1DHH / T1 - (T1DH / T1) * (T1DH / T1) + E0DHH / E0 -
                                    (E0DH / E0) * (E0DH / E0) -
                                    E1DHH / E1 + (E1DH / E1) * (E1DH / E1) - 1.0_rt / THETA2);
        Real EXP1TH = std::exp(-1.0_rt / THETA);
        Real C = (0.872496_rt + 0.025248_rt * EXP1TH) * E;
        Real CDH = 0.025248_rt * EXP1TH / THETA2 * E + C * EDH / E;
        Real CDHH = 0.025248_rt * EXP1TH / THETA2 * (EDH + (1.0_rt - 2.0_rt * THETA) / THETA2 * E) +
                    CDH * EDH / E + C * EDHH / E - C * (EDH / E) * (EDH / E);
        Real DISCR = std::sqrt(4.0_rt * E - D * D);
        Real DIDH = 0.5_rt / DISCR * (4.0_rt * EDH - 2.0_rt * D * DDH);
        Real DIDHH = (-std::pow((2.0_rt * EDH - D * DDH) / DISCR, 2) + 2.0_rt * EDHH -
                      DDH * DDH - D * DDHH) / DISCR;
        Real S1 = -C / E * GAME;
        Real S1H = CDH / C - EDH / E;
        Real S1DH = S1 * S1H;
        Real S1DHH = S1DH * S1H + S1 * (CDHH / C - (CDH / C) * (CDH / C) -
                                        EDHH / E + (EDH / E) * (EDH / E));
        Real S1DG = -C / E; // = > S1DGG = 0
        Real S1DHG = S1DG * (CDH / C - EDH / E);
        Real B2 = B - C * D / E;
        Real B2DH = BDH - (CDH * D + C * DDH) / E + C * D * EDH / (E * E);
        Real B2DHH = BDHH - (CDHH * D + 2.0_rt * CDH * DDH + C * DDHH) / E +
                     (2.0_rt * (CDH * D + C * DDH - C * D * EDH / E) * EDH +
                      C * D * EDHH) / (E * E);
        Real SQGE = std::sqrt(GAME);
        Real S2 = -2.0_rt / E * B2 * SQGE;
        Real S2H = B2DH / B2 - EDH / E;
        Real S2DH = S2 * S2H;
        Real S2DHH = S2DH * S2H + S2 * (B2DHH / B2 - (B2DH / B2) * (B2DH / B2) -
                                        EDHH / E + (EDH / E) * (EDH / E));
        Real S2DG = 0.5_rt * S2 / GAME;
        Real S2DGG = -0.5_rt * S2DG / GAME;
        Real S2DHG = 0.5_rt * S2DH / GAME;
        Real R3 = E * GAME + D * SQGE + 1.0_rt;
        Real R3DH = EDH * GAME + DDH * SQGE;
        Real R3DHH = EDHH * GAME + DDHH * SQGE;
        Real R3DG = E + 0.5_rt * D / SQGE;
        Real R3DGG = -0.25_rt * D / (GAME * SQGE);
        Real R3DHG = EDH + 0.5_rt * DDH / SQGE;
        Real B3 = A - C / E;
        Real B3DH = ADH - CDH / E + C * EDH / (E * E);
        Real B3DHH = ADHH - CDHH / E + (2.0_rt * CDH * EDH + C * EDHH) / (E * E) -
                     2.0_rt * C * EDH * EDH / (E * E * E);
        Real C3 = (D / E * B2 - B3) / E; // = D * B2 / (E * E) - B3 / E;
        Real C3DH = (DDH * B2 + D * B2DH + B3 * EDH) / (E * E) -
                    2.0_rt * D * B2 * EDH / (E * E * E) - B3DH / E;
        Real C3DHH = (-B3DHH +
                      (DDHH * B2 + 2.0_rt * DDH * B2DH + D * B2DHH +
                       B3DH * EDH + B3 * EDHH + B3DH * EDH) / E -
                      2.0_rt * ((DDH * B2 + D * B2DH + B3 * EDH + DDH * B2 + D * B2DH) * EDH +
                                D * B2 * EDHH) / (E * E) +
                      6.0_rt * D * B2 * EDH * EDH / (E * E * E)) / E;
        Real S3 = C3 * std::log(R3);
        Real S3DH = S3 * C3DH / C3 + C3 * R3DH / R3;
        Real S3DHH = (S3DH * C3DH + S3 * C3DHH) / C3 - S3 * (C3DH / C3) * (C3DH / C3) +
                     (C3DH * R3DH + C3 * R3DHH) / R3 - C3 * (R3DH / R3) * (R3DH / R3);
        Real S3DG = C3 * R3DG / R3;
        Real S3DGG = C3 * (R3DGG / R3 - (R3DG / R3) * (R3DG / R3));
        Real S3DHG = (C3DH * R3DG + C3 * R3DHG) / R3 - C3 * R3DG * R3DH / (R3 * R3);
        Real B4 = 2.0_rt - D * D / E;
        Real B4DH = EDH * (D / E) * (D / E) - 2.0_rt * D * DDH / E;
        Real B4DHH = EDHH * (D / E) * (D / E) + 2.0_rt * EDH * (D / E) * (D / E) * (DDH / D - EDH / E) -
                     2.0_rt * (DDH * DDH + D * DDHH) / E + 2.0_rt * D * DDH * EDH / (E * E);
        Real C4 = 2.0_rt * E * SQGE + D;
        Real C4DH = 2.0_rt * EDH * SQGE + DDH;
        Real C4DHH = 2.0_rt * EDHH * SQGE + DDHH;
        Real C4DG = E / SQGE;
        Real C4DGG = -0.5_rt * E / (GAME * SQGE);
        Real C4DHG = EDH / SQGE;
        Real S4A = 2.0_rt / E / DISCR;
        Real S4AH = EDH / E + DIDH / DISCR;
        Real S4ADH = -S4A * S4AH;
        Real S4ADHH = -S4ADH * S4AH -
                       S4A * (EDHH / E - (EDH / E) * (EDH / E) + DIDHH / DISCR -
                              (DIDH / DISCR) * (DIDH / DISCR));
        Real S4B = D * B3 + B4 * B2;
        Real S4BDH = DDH * B3 + D * B3DH + B4DH * B2 + B4 * B2DH;
        Real S4BDHH = DDHH * B3 + 2.0_rt * DDH * B3DH + D * B3DHH + B4DHH * B2 +
                      2.0_rt * B4DH * B2DH + B4 * B2DHH;
        Real S4C = std::atan(C4 / DISCR) - std::atan(D / DISCR);
        Real UP1 = C4DH * DISCR - C4 * DIDH;
        Real DN1 = DISCR * DISCR + C4 * C4;
        Real UP2 = DDH * DISCR - D * DIDH;
        Real DN2 = DISCR * DISCR + D * D;
        Real S4CDH = UP1 / DN1 - UP2 / DN2;
        Real S4CDHH = (C4DHH * DISCR - C4 * DIDHH) / DN1 -
                      UP1 * 2.0_rt * (DISCR * DIDH + C4 * C4DH) / (DN1 * DN1) -
                      (DDHH * DISCR - D * DIDHH) / DN2 + UP2 * 2.0_rt *
                      (DISCR * DIDH + D * DDH) / (DN2 * DN2);
        Real S4CDG = C4DG * DISCR / DN1;
        Real S4CDGG = C4DGG * DISCR / DN1 - 2.0_rt * C4 * DISCR * (C4DG / DN1) * (C4DG / DN1);
        Real S4CDHG = (C4DHG * DISCR + C4DG * DIDH -
                       C4DG * DISCR / DN1 * 2.0_rt * (DISCR * DIDH + C4 * C4DH)) / DN1;
        Real S4 = S4A * S4B * S4C;
        Real S4DH = S4ADH * S4B * S4C + S4A * S4BDH * S4C + S4A * S4B * S4CDH;
        Real S4DHH = S4ADHH * S4B * S4C + S4A * S4BDHH * S4C + S4A * S4B * S4CDHH +
                     2.0_rt * (S4ADH * S4BDH * S4C + S4ADH * S4B * S4CDH + S4A * S4BDH * S4CDH);
        Real S4DG = S4A * S4B * S4CDG;
        Real S4DGG = S4A * S4B * S4CDGG;
        Real S4DHG = S4A * S4B * S4CDHG + S4CDG * (S4ADH * S4B + S4A * S4BDH);

        FXC = S1 + S2 + S3 + S4;
        Real FXCDH = S1DH + S2DH + S3DH + S4DH;
        Real FXCDG = S1DG + S2DG + S3DG + S4DG;
        Real FXCDHH = S1DHH + S2DHH + S3DHH + S4DHH;
        Real FXCDGG = S2DGG + S3DGG + S4DGG;
        Real FXCDHG = S1DHG + S2DHG + S3DHG + S4DHG;
        PXC = (GAME * FXCDG - 2.0_rt * THETA * FXCDH) / 3.0_rt;
        UXC = GAME * FXCDG - THETA * FXCDH;
        SXC = (GAME * S2DG - S2 + GAME * S3DG - S3 + S4A * S4B * (GAME * S4CDG - S4C)) -
              THETA * FXCDH;
        if (std::abs(SXC) < EPS * std::abs(THETA * FXCDH)) {
            SXC = 0.0_rt; // accuracy loss
        }
        CVXC = 2.0_rt * THETA * (GAME * FXCDHG - FXCDH) - THETA * THETA * FXCDHH - GAME * GAME * FXCDGG;
        if (std::abs(CVXC) < EPS * std::abs(GAME * GAME * FXCDGG)) {
            CVXC = 0.0_rt; // accuracy
        }
        Real PDLH = THETA * (GAME * FXCDHG - 2.0_rt * FXCDH - 2.0_rt * THETA * FXCDHH) / 3.0_rt;
        Real PDLG = GAME * (FXCDG + GAME * FXCDGG - 2.0_rt * THETA * FXCDHG) / 3.0_rt;
        PDRXC = PXC + (PDLG - 2.0_rt * PDLH) / 3.0_rt;
        PDTXC = GAME * (THETA * FXCDHG - GAME * FXCDGG / 3.0_rt) -
                THETA * (FXCDH / 0.75_rt + THETA * FXCDHH / 1.5_rt);
    }

    void subfermj (Real CMU1,
                   Real& CJ00, Real& CJ10, Real& CJ20,
                   Real& CJ01, Real& CJ11, Real& CJ21,
                   Real& CJ02, Real& CJ12, Real& CJ22,
                   Real& CJ03, Real& CJ13, Real& CJ23,
                   Real& CJ04, Real& CJ14, Real& CJ24, Real& CJ05)
    {
        // Version 17.11.11
        // corrected 04.03.21
        // Supplement to SOMMERF
        const Real EPS = 1.e-4_rt; // inserted 04.03.21
        if (CMU1 <= 0.0_rt) {
            printf("SUBFERMJ: small CHI\n");
            exit(1);
        }

        Real CMU = 1.0_rt + CMU1;
        Real X0 = std::sqrt(CMU1 * (2.0_rt + CMU1));
        Real X3 = X0 * X0 * X0;
        Real X5 = X3 * X0 * X0;
        Real X7 = X5 * X0 * X0;
        if (X0 < EPS) {
            CJ00 = X3 / 3.0_rt;
            CJ10 = 0.1_rt * X5;
            CJ20 = X7 / 28.0_rt;
        }
        else {
            Real CL = std::log(X0 + CMU);
            CJ00 = 0.5_rt * (X0 * CMU - CL); // J_{1/2}^0
            CJ10 = X3 / 3.0_rt - CJ00; // J_{3/2}^0
            CJ20 = (0.75_rt * CMU - 2.0_rt) / 3.0_rt * X3 + 1.25_rt * CJ00; // J_{5/2}^0
        }

        CJ01 = X0; // J_{1/2}^1
        CJ11 = CJ01 * CMU1; // J_{3/2}^1
        CJ21 = CJ11 * CMU1; // J_{5/2}^1
        Real RCJ02 = CMU / X0; // J_{1/2}^2
        CJ12 = CMU1 / X0 * (3.0_rt + 2.0_rt * CMU1); // J_{3/2}^2
        CJ22 = CMU1 * CMU1 / X0 * (5.0_rt + 3.0_rt * CMU1); // J_{5/2}^2
        CJ03 = -1.0_rt / X3; // J_{1/2}^3
        CJ13 = CMU1 / X3 * (2.0_rt * CMU1 * CMU1 + 6.0_rt * CMU1 + 3.0_rt);
        CJ23 = CMU1 * CMU1 / X3 * (6.0_rt * CMU1 * CMU1 + 2.0e1_rt * CMU1 + 1.5e1_rt);
        CJ04 = 3.0_rt * CMU / X5;
        CJ14 = -3.0_rt * CMU1 / X5;
        CJ24 = CMU1 * CMU1 / X5 * (6.0_rt * CMU1 * CMU1 * CMU1 + 3.0e1_rt * CMU1 * CMU1 +
                                   45.0_rt * CMU1 + 15.0_rt);
        CJ05 = (-12.0_rt * CMU1 * CMU1 - 24.0_rt * CMU1 - 15.0_rt) / (X7);
    }

    void sommerf (Real TEMR, Real CHI,
                  Real& W0, Real& W0DX, Real& W0DT, Real& W0DXX, Real& W0DTT, Real& W0DXT,
                  Real& W1, Real& W1DX, Real& W1DT, Real& W1DXX, Real& W1DTT, Real& W1DXT,
                  Real& W2, Real& W2DX, Real& W2DT, Real& W2DXX, Real& W2DTT, Real& W2DXT,
                  Real& W0XXX, Real& W0XTT, Real& W0XXT)
    {
        // Version 17.11.11
        // Sommerfeld expansion for the Fermi-Dirac integrals
        // Input: TEMR=T/mc^2; CHI=(\mu-mc^2)/T
        // Output: Wk - Fermi-Dirac integral of the order k+1/2
        //         WkDX=dWk/dCHI, WkDT = dWk/dT, WkDXX=d^2 Wk / d CHI^2,
        //         WkDTT=d^2 Wk / d T^2, WkDXT=d^2 Wk /dCHIdT,
        //         W0XXX=d^3 W0 / d CHI^3, W0XTT=d^3 W0 /(d CHI d^2 T),
        //         W0XXT=d^3 W0 /dCHI^2 dT
        // [Draft source: yellow book pages 124-127]

        const Real PI = 3.141592653_rt;
        const Real PI2 = PI * PI;

        if (CHI < 0.5_rt) {
            printf("SOMMERF: non-degenerate (small CHI)\n");
            exit(1);
        }

        if (TEMR <= 0.0_rt) {
            printf("SOMMERF: T < 0\n");
            exit(1);
        }

        Real CMU1 = CHI * TEMR; // chemical potential in rel.units
        Real CMU = 1.0_rt + CMU1;

        Real CJ00, CJ10, CJ20;
        Real CJ01, CJ11, CJ21;
        Real CJ02, CJ12, CJ22;
        Real CJ03, CJ13, CJ23;
        Real CJ04, CJ14, CJ24;
        Real CJ05;

        subfermj(CMU1,
                 CJ00, CJ10, CJ20,
                 CJ01, CJ11, CJ21,
                 CJ02, CJ12, CJ22,
                 CJ03, CJ13, CJ23,
                 CJ04, CJ14, CJ24, CJ05);

        Real PIT26 = (PI * TEMR)*(PI * TEMR) / 6.0_rt;
        Real CN0 = std::sqrt(0.5_rt / TEMR) / TEMR;
        Real CN1 = CN0 / TEMR;
        Real CN2 = CN1 / TEMR;
        W0 = CN0 * (CJ00 + PIT26 * CJ02); // + CN0 * PITAU4 * CJ04
        W1 = CN1 * (CJ10 + PIT26 * CJ12); // + CN1 * PITAU4 * CJ14
        W2 = CN2 * (CJ20 + PIT26 * CJ22); // + CN2 * PITAU4 * CJ24
        W0DX = CN0 * TEMR * (CJ01 + PIT26 * CJ03); //  + CN0 * PITAU4 * CJ05
        W1DX = CN0 * (CJ11 + PIT26 * CJ13);
        W2DX = CN1 * (CJ21 + PIT26 * CJ23);
        W0DT = CN1 * (CMU1 * CJ01 - 1.5_rt * CJ00 + PIT26 * (CMU1 * CJ03 + 0.5_rt * CJ02));
        W1DT = CN2 * (CMU1 * CJ11 - 2.5_rt * CJ10 + PIT26 * (CMU1 * CJ13 - 0.5_rt * CJ12));
        W2DT = CN2 / TEMR * (CMU1 * CJ21 - 3.5_rt * CJ20 + PIT26 * (CMU1 * CJ23 - 1.5_rt * CJ22));
        W0DXX = CN0 * TEMR * TEMR * (CJ02 + PIT26 * CJ04);
        W1DXX = CN0 * TEMR * (CJ12 + PIT26 * CJ14);
        W2DXX = CN0 * (CJ22 + PIT26 * CJ24);
        W0DXT = CN0 * (CMU1 * CJ02 - 0.5_rt * CJ01 + PIT26 * (CMU1 * CJ04 + 1.5_rt * CJ03));
        W1DXT = CN1 * (CMU1 * CJ12 - 1.5_rt * CJ11 + PIT26 * (CMU1 * CJ14 + 0.5_rt * CJ13));
        W2DXT = CN2 * (CMU1 * CJ22 - 2.5_rt * CJ21 + PIT26 * (CMU1 * CJ24 - 0.5_rt * CJ23));
        W0DTT = CN2 * (3.75_rt * CJ00 - 3.0_rt * CMU1 * CJ01 + CMU1 * CMU1 * CJ02 +
                       PIT26 * (-0.25_rt * CJ02 + CMU1 * CJ03 + CMU1 * CMU1 * CJ04));
        W1DTT = CN2 / TEMR * (8.75_rt * CJ10 - 5.0_rt * CMU1 * CJ11 + CMU1 * CMU1 * CJ12 +
                              PIT26 * (0.75_rt * CJ12 - CMU1 * CJ13 + CMU1 * CMU1 * CJ14));
        W2DTT = CN2 / TEMR * TEMR * (15.75_rt * CJ20 - 7.0_rt * CMU1 * CJ21 + CMU1 * CMU1 * CJ22 +
                                     PIT26 * (3.75_rt * CJ22 - 3.0_rt * CMU1 * CJ23 + CMU1 * CMU1 * CJ24));
        W0XXX = CN0 * TEMR * TEMR * TEMR * (CJ03 + PIT26 * CJ05);
        W0XXT = CN0 * TEMR * (CMU1 * CJ03 + 0.5_rt * CJ02 + PIT26 * (CMU1 * CJ05 + 2.5_rt * CJ04));
        W0XTT = CN1 * (0.75_rt * CJ01 - CMU1 * CJ02 + CMU1 * CMU1 * CJ03 +
                       PIT26 * (0.75_rt * CJ03 + 3.0_rt * CMU1 * CJ04 + CMU1 * CMU1 * CJ05));
    }

    void elect11b(Real TEMP, Real CHI,
                  Real& DENS, Real& FEid, Real& PEid, Real& UEid,
                  Real& SEid, Real& CVE, Real& CHITE, Real& CHIRE,
                  Real& DlnDH, Real& DlnDT, Real& DlnDHH,
                  Real& DlnDTT, Real& DlnDHT)
    {
        // Version 17.11.11
        // Stems from ELECT9b v.19.01.10, Diff. - additional output.
        // Sommerfeld expansion at very large CHI.

        const Real BOHR = 137.036_rt;
        const Real PI = 3.141592653_rt;
        const Real PI2 = PI * PI;
        const Real BOHR2 = BOHR * BOHR;
        const Real BOHR3 = BOHR2 * BOHR; // cleaned 15/6

        Real TEMR = TEMP / BOHR2; // T in rel.units ( = T/mc^2)
        Real EF = CHI * TEMR; // Fermi energy in mc^2 - zeroth aprox.  =  CMU1
        Real DeltaEF = PI2 * TEMR * TEMR / 6.0_rt * (1.0_rt + 2.0_rt * EF * (2.0_rt + EF)) /
                       (EF * (1.0_rt + EF) * (2.0_rt + EF)); // corr. [p.125, equiv.Eq.(6) of PC'10]
        EF = EF + DeltaEF; // corrected Fermi energy (14.02.09)
        Real G = 1.0_rt + EF; // electron Lorentz-factor

        Real PF, F, DF, P, DP;

        if (EF > 1.e-5_rt) { // relativistic expansion (Yak.&Shal.'89)
            PF = std::sqrt(G * G - 1.0_rt); // Fermi momentum [rel.un. = mc]
            F = (PF * (1.0_rt + 2.0_rt * PF * PF) * G - PF * PF * PF / .375_rt - std::log(PF + G)) / 8.0_rt / PI2; // F/V
            DF = -TEMR * TEMR * PF * G / 6.0_rt; // thermal correction to F/V
            P = (PF * G * (PF * PF / 1.5_rt - 1.0_rt) + std::log(PF + G)) / 8.0_rt / PI2; // P(T = 0)
            DP = TEMR * TEMR * PF * (PF * PF + 2.0_rt) / G / 18.0_rt; // thermal correction to P
            CVE = PI2 * TEMR * G / (PF * PF);
        }
        else { // nonrelativistic limit
            PF = std::sqrt(2.0_rt * EF);
            F = (PF * PF * PF * PF * PF) * 0.1_rt / PI2;
            DF = -TEMR * TEMR * PF / 6.0_rt;
            P = F / 1.5_rt;
            DP = TEMR * TEMR * PF / 9.0_rt;
            CVE = PI2 * TEMR / EF / 2.0_rt;
        }

        F = F + DF;
        P = P + DP;
        Real S = -2.0_rt * DF; // entropy per unit volume [rel.un.]
        Real U = F + S;
        CHIRE = (PF * PF * PF * PF * PF) / (9.0_rt * PI2 * P * G);
        CHITE = 2.0_rt * DP / P;
        Real DENR = PF * PF * PF / 3.0_rt / PI2; // n_e [rel.un. = \Compton^{-3}]
        DENS = DENR * BOHR3; // conversion to a.u.( = \Bohr_radius^{-3})

        // derivatives over chi at constant T and T at constant chi:
        Real TPI = TEMR * std::sqrt(2.0_rt * TEMR) / PI2; // common pre-factor

        Real W0, W0DX, W0DT, W0DXX, W0DTT, W0DXT;
        Real W1, W1DX, W1DT, W1DXX, W1DTT, W1DXT;
        Real W2, W2DX, W2DT, W2DXX, W2DTT, W2DXT;
        Real W0XXX, W0XTT, W0XXT;

        sommerf(TEMR, CHI,
                W0, W0DX, W0DT, W0DXX, W0DTT, W0DXT,
                W1, W1DX, W1DT, W1DXX, W1DTT, W1DXT,
                W2, W2DX, W2DT, W2DXX, W2DTT, W2DXT,
                W0XXX, W0XTT, W0XXT);

        Real dndH = TPI * (W0DX + TEMR * W1DX); // (d n_e/d\chi)_T
        Real dndT = TPI * (1.5_rt * W0 / TEMR + 2.5 * W1 + W0DT + TEMR * W1DT); // (d n_e/dT)_\chi
        Real dndHH = TPI * (W0DXX + TEMR * W1DXX); // (d^2 n_e/d\chi)_T
        Real dndTT = TPI * (0.75_rt * W0 / TEMR * TEMR + 3. * W0DT / TEMR + W0DTT +
                            3.75 * W1 / TEMR + 5. * W1DT + TEMR * W1DTT);
        Real dndHT = TPI * (1.5_rt * W0DX / TEMR + W0DXT + 2.5 * W1DX + TEMR * W1DXT);

        DlnDH = dndH / DENR; // (d ln n_e/d\chi)_T
        DlnDT = dndT * TEMR / DENR; // (d ln n_e/d ln T)_\chi
        DlnDHH = dndHH / DENR - DlnDH * DlnDH; // (d^2 ln n_e/d\chi^2)_T
        DlnDTT = TEMR * TEMR / DENR * dndTT + DlnDT - DlnDT * DlnDT; // d^2 ln n_e/d ln T^2
        DlnDHT = TEMR / DENR * (dndHT - dndT * DlnDH); // d^2 ln n_e/d\chi d ln T

        Real DT = DENR * TEMR;
        PEid = P / DT;
        UEid = U / DT;
        FEid = F / DT;
        SEid = S / DT;

        // Empirical corrections of 16.02.09:
        Real D1 = DeltaEF / EF;
        Real D2 = D1 * (4.0_rt - 2.0_rt * (PF / G));
        CVE = CVE / (1.0_rt + D2);
        SEid = SEid / (1.0_rt + D1);
        CHITE = CHITE / (1.0_rt + D2);
    }

    void elect11a(Real TEMP, Real CHI,
                  Real& DENS, Real& FEid, Real& PEid, Real& UEid,
                  Real& SEid, Real& CVE, Real& CHITE, Real& CHIRE,
                  Real& DlnDH, Real& DlnDT, Real& DlnDHH, Real& DlnDTT,
                  Real& DlnDHT)
    {
        // Version 16.11.11
        // This is THE FIRST PART of ELECT9 v.04.03.09.
        const Real BOHR = 137.036_rt;
        const Real PI = 3.141592653_rt;
        const Real PI2 = PI * PI;
        const Real BOHR2 = BOHR * BOHR;
        const Real BOHR3 = BOHR2 * BOHR; // cleaned 15/6

        Real TEMR = TEMP / BOHR2; // T in rel.units (=T/mc^2)

        Real W0, W0DX, W0DT, W0DXX, W0DTT, W0DXT;
        Real W1, W1DX, W1DT, W1DXX, W1DTT, W1DXT;
        Real W2, W2DX, W2DT, W2DXX, W2DTT, W2DXT;
        Real W0XXX, W0XTT, W0XXT;

        blin9(TEMR, CHI,
              W0, W0DX, W0DT, W0DXX, W0DTT, W0DXT,
              W1, W1DX, W1DT, W1DXX, W1DTT, W1DXT,
              W2, W2DX, W2DT, W2DXX, W2DTT, W2DXT,
              W0XXX, W0XTT, W0XXT);

        Real TPI = TEMR * std::sqrt(2.0_rt * TEMR) / PI2; // common pre-factor
        Real DENR = TPI * (W1 * TEMR + W0);
        Real PR = TEMR * TPI / 3.0_rt * (W2 * TEMR + 2.0_rt * W1);
        Real U = TEMR * TPI * (W2 * TEMR + W1);

        // (these are density, pressure, and internal energy in the rel.units)
        PEid = PR / (DENR * TEMR);
        UEid = U / (DENR * TEMR);
        FEid = CHI - PEid;
        DENS = DENR * BOHR3; // converts from rel.units to a.u.
        SEid = UEid - FEid;

        // derivatives over T at constant chi:
        Real dndT = TPI * (1.5_rt * W0 / TEMR + 2.5_rt * W1 + W0DT + TEMR * W1DT); // (d n_e/dT)_\chi
        Real dPdT = TPI / 3.0_rt * (5.0_rt * W1 + 2.0_rt * TEMR * W1DT + 3.5_rt * TEMR * W2 + TEMR * TEMR * W2DT); //dP/dT
        Real dUdT = TPI * (2.5_rt * W1 + TEMR * W1DT + 3.5_rt * TEMR * W2 + TEMR * TEMR * W2DT); //dU/dT_\chi

        // derivatives over chi at constant T and second derivatives:
        Real dndH = TPI * (W0DX + TEMR * W1DX); // (d n_e/d\chi)_T
        Real dndHH = TPI * (W0DXX + TEMR * W1DXX); // (d^2 n_e/d\chi)_T
        Real dndTT = TPI * (0.75_rt * W0 / TEMR * TEMR + 3.0_rt * W0DT / TEMR + W0DTT +
                            3.75_rt * W1 / TEMR + 5.0_rt * W1DT + TEMR * W1DTT);
        Real dndHT = TPI * (1.5_rt * W0DX / TEMR + W0DXT + 2.5_rt * W1DX + TEMR * W1DXT);

        DlnDH = dndH / DENR; // (d ln n_e/d\chi)_T
        DlnDT = dndT * TEMR / DENR; // (d ln n_e/d ln T)_\chi
        DlnDHH = dndHH / DENR - DlnDH * DlnDH; // (d^2 ln n_e/d\chi^2)_T
        DlnDTT = TEMR * TEMR / DENR * dndTT + DlnDT - DlnDT * DlnDT; // d^2 ln n_e/d ln T^2
        DlnDHT = TEMR / DENR * (dndHT - dndT * DlnDH); // d^2 ln n_e/d\chi d ln T
        Real dPdH = TPI / 3.0_rt * TEMR * (2.0_rt * W1DX + TEMR * W2DX); // (d P_e/d\chi)_T
        Real dUdH = TPI * TEMR * (W1DX + TEMR * W2DX); // (d U_e/d\chi)_T
        CVE = (dUdT - dUdH * dndT / dndH) / DENR;
        CHITE = TEMR / PR * (dPdT - dPdH * dndT / dndH);
        CHIRE = DENR / PR * dPdH / dndH; // (dndH * TEMR * PEid) // DENS / PRE * dPdH / dndH
    }

    void elect11 (double TEMP, double CHI,
                  Real& DENS, Real& FEid, Real& PEid, Real& UEid,
                  Real& SEid, Real& CVE, Real& CHITE, Real& CHIRE,
                  Real& DlnDH, Real& DlnDT, Real& DlnDHH, Real& DlnDTT,
                  Real& DlnDHT)
    {
        // Version 17.11.11
        // safeguard against huge (-CHI) values is added 27.05.17
        // ELECT9 v.04.03.09 + smooth match of two fits at chi >> 1 + add.outputs
        // Compared to ELECTRON v.06.07.00, this S/R is completely rewritten:
        //        numerical differentiation is avoided now.
        // Compared to ELECT7 v.06.06.07,
        //    - call BLIN7 is changed to call BLIN9,
        //    - Sommerfeld expansion is used at chi >~ 28 i.o. 1.e4
        //    - Sommerfeld expansion is corrected: introduced DeltaEF, D1 and D2.
        // Ideal electron-gas EOS.
        // Input: TEMP - T [a.u.], CHI=\mu/T
        // Output: DENS - electron number density n_e [a.u.],
        //         FEid - free energy / N_e kT, UEid - internal energy / N_e kT,
        //         PEid - pressure (P_e) / n_e kT, SEid - entropy / N_e k,
        //         CVE - heat capacity / N_e k,
        //         CHITE=(d ln P_e/d ln T)_V, CHIRE=(d ln P_e/d ln n_e)_T
        //         DlnDH=(d ln n_e/d CHI)_T = (T/n_e) (d n_e/d\mu)_T
        //         DlnDT=(d ln n_e/d ln T)_CHI
        //         DlnDHH=(d^2 ln n_e/d CHI^2)_T
        //         DlnDTT=(d^2 ln n_e/d (ln T)^2)_CHI
        //         DlnDHT=d^2 ln n_e/d (ln T) d CHI

        const Real CHI2 = 28.0_rt;
        const Real XMAX = 20.0_rt;
        const Real DCHI2 = CHI2 - 1.0_rt;
        const Real XSCAL2 = XMAX / DCHI2;

        if (CHI < -1.e2_rt) {
            printf("ELECT11: too large negative CHI\n"); // 27.05.17
            exit(1);
        }

        Real X2 = (CHI - CHI2) * XSCAL2;
        if (X2 < -XMAX) {
            elect11a(TEMP, CHI,
                     DENS, FEid, PEid, UEid, SEid, CVE, CHITE, CHIRE,
                     DlnDH, DlnDT, DlnDHH, DlnDTT, DlnDHT);
        }
        else if (X2 > XMAX) {
            elect11b(TEMP, CHI,
                     DENS, FEid, PEid, UEid, SEid, CVE, CHITE, CHIRE,
                     DlnDH, DlnDT, DlnDHH, DlnDTT, DlnDHT);
        }
        else {
            Real FP, FM;
            fermi10(X2, XMAX, FP, FM);

            Real DENSa, FEida, PEida, UEida, SEida, CVEa, CHITEa, CHIREa;
            Real DlnDHa, DlnDTa, DlnDHHa, DlnDTTa, DlnDHTa;

            elect11a(TEMP, CHI,
                     DENSa, FEida, PEida, UEida, SEida, CVEa, CHITEa, CHIREa,
                     DlnDHa, DlnDTa, DlnDHHa, DlnDTTa, DlnDHTa);

            Real DENSb, FEidb, PEidb, UEidb, SEidb, CVEb, CHITEb, CHIREb;
            Real DlnDHb, DlnDTb, DlnDHHb, DlnDTTb, DlnDHTb;

            elect11b(TEMP, CHI,
                     DENSb, FEidb, PEidb, UEidb, SEidb, CVEb, CHITEb, CHIREb,
                     DlnDHb, DlnDTb, DlnDHHb, DlnDTTb, DlnDHTb);

            DENS = DENSa * FP + DENSb * FM;
            FEid = FEida * FP + FEidb * FM;
            PEid = PEida * FP + PEidb * FM;
            UEid = UEida * FP + UEidb * FM;
            SEid = SEida * FP + SEidb * FM;
            CVE = CVEa * FP + CVEb * FM;
            CHITE = CHITEa * FP + CHITEb * FM;
            CHIRE = CHIREa * FP + CHIREb * FM;
            DlnDH = DlnDHa * FP + DlnDHb * FM;
            DlnDT = DlnDTa * FP + DlnDTb * FM;
            DlnDHH = DlnDHHa * FP + DlnDHHb * FM;
            DlnDHT = DlnDHTa * FP + DlnDHTb * FM;
            DlnDTT = DlnDTTa * FP + DlnDTTb * FM;
        }
    }

    void fscrsol8 (Real RS, Real GAMI, Real ZNUCL, Real TPT,
                   Real& FSCR, Real& USCR, Real& PSCR, Real& S_SCR,
                   Real& CVSCR, Real& PDTSCR, Real& PDRSCR)
    {
        // Version 28.05.08
        // undefined zero variable Q1DXG is wiped out 21.06.10
        // accuracy-loss safeguard added 10.08.16
        // safequard against Zion < 1 added 27.05.17
        // Fit to the el.-ion screening in bcc or fcc Coulomb solid
        // Stems from FSCRsol8 v.09.06.07. Included a check for RS = 0.
        //   INPUT: RS - el. density parameter, GAMI - ion coupling parameter,
        //          ZNUCL - ion charge, TPT = T_p/T - ion quantum parameter
        //   OUTPUT: FSCR - screening (e-i) free energy per kT per 1 ion,
        //           USCR - internal energy per kT per 1 ion (screen.contrib.)
        //           PSCR - pressure divided by (n_i kT) (screen.contrib.)
        //           S_SCR - screening entropy contribution / (N_i k)
        //           CVSCR - heat capacity per 1 ion (screen.contrib.)
        //           PDTSCR,PDRSCR  =  PSCR + d PSCR / d ln(T,\rho)

        const Real C13 = 1.0_rt / 3.0_rt;
        const Real ENAT = 2.7182818285_rt;
        const Real TINY = 1.e-19_rt;

        const Real AP[4] = {1.1866_rt, 0.684_rt, 17.9_rt, 41.5_rt};
        const Real PX = 0.205_rt; // for bcc lattice

        if (RS < 0.0_rt) {
            printf("FSCRliq8: RS < 0\n");
            exit(1);
        }

        if (RS < TINY) {
            FSCR = 0.0_rt;
            USCR = 0.0_rt;
            PSCR = 0.0_rt;
            S_SCR = 0.0_rt;
            CVSCR = 0.0_rt;
            PDTSCR = 0.0_rt;
            PDRSCR = 0.0_rt;
            return;
        }

        Real Zion = ZNUCL;
        if (Zion < 1.0_rt) { // 27.05.17
            Zion = 1.0_rt;
        }

        Real XSR = 0.0140047_rt / RS; // relativity parameter
        Real Z13 = std::pow(Zion, C13);
        Real P1 = 0.00352_rt * (1.0_rt - AP[0] / std::pow(Zion, 0.267_rt) + 0.27_rt / Zion);
        Real P2 = 1.0_rt + 2.25_rt / Z13 *
                  (1.0_rt + AP[1] * (Zion * Zion * Zion * Zion * Zion) +
                   0.222_rt * (Zion * Zion * Zion * Zion * Zion * Zion)) /
                  (1.0_rt + .222 * Zion * Zion * Zion * Zion * Zion * Zion);
        Real ZLN = std::log(Zion);
        Real Finf = std::sqrt(P2 / (XSR * XSR) + 1.0_rt) * Z13 * Z13 * P1; // The TF limit
        Real FinfX = -P2 / ((P2 + XSR * XSR) * XSR);
        Real FinfDX = Finf * FinfX;
        Real FinfDXX = FinfDX * FinfX - FinfDX * (P2 + 3.0_rt * XSR * XSR) / ((P2 + XSR * XSR) * XSR);
        Real R1 = AP[3] / (1.0_rt + ZLN);
        Real R2 = 0.395_rt * ZLN + .347 / Zion / std::sqrt(Zion);
        Real R3 = 1.0_rt / (1.0_rt + ZLN * std::sqrt(ZLN) * 0.01_rt + 0.097_rt / (Zion * Zion));
        Real Q1U = R1 + AP[2] * XSR * XSR;
        Real Q1D = 1.0_rt + R2 * XSR * XSR;
        Real Q1 = Q1U / Q1D;
        Real Q1X = 2.0_rt * XSR * (AP[2] / Q1U - R2 / Q1D);
        Real Q1XDX = Q1X / XSR + 4.0_rt * XSR * XSR * ((R2 / Q1D) * (R2 / Q1D) - (AP[2] / Q1U) * (AP[2] / Q1U));
        Real Q1DX = Q1 * Q1X;
        Real Q1DXX = Q1DX * Q1X + Q1 * Q1XDX;

        Real SUP, SUPDX, SUPDG, SUPDXX, SUPDGG, SUPDXG;

        // New quantum factor, in order to suppress CVSCR at TPT >> 1
        if (TPT < 6.0_rt / PX) {
            Real Y0 = (PX * TPT) * (PX * TPT);
            Real Y0DX = Y0 / XSR;
            Real Y0DG = 2.0_rt * Y0 / GAMI;
            Real Y0DXX = 0.0_rt;
            Real Y0DGG = Y0DG / GAMI;
            Real Y0DXG = Y0DG / XSR;
            Real Y1 = std::exp(Y0);
            Real Y1DX = Y1 * Y0DX;
            Real Y1DG = Y1 * Y0DG;
            Real Y1DXX = Y1 * (Y0DX * Y0DX + Y0DXX);
            Real Y1DGG = Y1 * (Y0DG * Y0DG + Y0DGG);
            Real Y1DXG = Y1 * (Y0DX * Y0DG + Y0DXG);
            Real SA = 1.0_rt + Y1;
            Real SUPA = std::log(SA);
            Real SUPADX = Y1DX / SA;
            Real SUPADG = Y1DG / SA;
            Real SUPADXX = (Y1DXX - Y1DX * Y1DX / SA) / SA;
            Real SUPADGG = (Y1DGG - Y1DG * Y1DG / SA) / SA;
            Real SUPADXG = (Y1DXG - Y1DX * Y1DG / SA) / SA;
            Real EM2 = ENAT - 2.0_rt;
            Real SB = ENAT - EM2 / Y1;
            Real SUPB = std::log(SB);
            Real EM2Y1 = EM2 / (Y1 * Y1 * SB);
            Real SUPBDX = EM2Y1 * Y1DX;
            Real SUPBDG = EM2Y1 * Y1DG;
            Real SUPBDXX = EM2Y1 * (Y1DXX - 2.0_rt * Y1DX * Y1DX / Y1 - Y1DX * SUPBDX);
            Real SUPBDGG = EM2Y1 * (Y1DGG - 2.0_rt * Y1DG * Y1DG / Y1 - Y1DG * SUPBDG);
            Real SUPBDXG = EM2Y1 * (Y1DXG - 2.0_rt * Y1DX * Y1DG / Y1 - Y1DG * SUPBDX);

            SUP = std::sqrt(SUPA / SUPB);
            Real SUPX = 0.5_rt * (SUPADX / SUPA - SUPBDX / SUPB);
            SUPDX = SUP * SUPX;
            Real SUPG = 0.5_rt * (SUPADG / SUPA - SUPBDG / SUPB);
            SUPDG = SUP * SUPG;
            SUPDXX = SUPDX * SUPX +
                     SUP * 0.5_rt * (SUPADXX / SUPA - (SUPADX / SUPA) * (SUPADX / SUPA) -
                                     SUPBDXX / SUPB + (SUPBDX / SUPB) * (SUPBDX / SUPB));
            SUPDGG = SUPDG * SUPG +
                     SUP * 0.5_rt * (SUPADGG / SUPA - (SUPADG / SUPA) * (SUPADG / SUPA) -
                                     SUPBDGG / SUPB + (SUPBDG / SUPB) * (SUPBDG / SUPB));
            SUPDXG = SUPDX * SUPG +
                     SUP * 0.5_rt * ((SUPADXG - SUPADX * SUPADG / SUPA) / SUPA -
                                     (SUPBDXG - SUPBDX * SUPBDG / SUPB) / SUPB);
        }
        else {
            SUP = PX * TPT;
            SUPDX = 0.5_rt * PX * TPT / XSR;
            SUPDG = PX * TPT / GAMI;
            SUPDXX =  - 0.5_rt * SUPDX / XSR;
            SUPDGG = 0.0_rt;
            SUPDXG = SUPDX / GAMI;
        }

        Real GR3 = std::pow(GAMI / SUP, R3);
        Real GR3X = -R3 * SUPDX / SUP;
        Real GR3DX = GR3 * GR3X;
        Real GR3DXX = GR3DX * GR3X - R3 * GR3 * (SUPDXX / SUP - (SUPDX / SUP) * (SUPDX / SUP));
        Real GR3G = R3 * (1.0_rt / GAMI - SUPDG / SUP);
        Real GR3DG = GR3 * GR3G;
        Real GR3DGG = GR3DG * GR3G + GR3 * R3 * ((SUPDG / SUP) * (SUPDG / SUP) - SUPDGG / SUP - 1.0_rt / (GAMI * GAMI));
        Real GR3DXG = GR3DG * GR3X + GR3 * R3 * (SUPDX * SUPDG / (SUP * SUP) - SUPDXG / SUP);
        Real W = 1.0_rt + Q1 / GR3;
        Real WDX = Q1DX / GR3 - Q1 * GR3DX / (GR3 * GR3);
        Real WDG = -Q1 * GR3DG / (GR3 * GR3);
        Real WDXX = Q1DXX / GR3 -
                    (2.0_rt * Q1DX * GR3DX + Q1 * (GR3DXX - 2.0_rt * GR3DX * GR3DX / GR3)) / (GR3 * GR3);
        Real WDGG = Q1 * (2.0_rt * GR3DG * GR3DG / GR3 - GR3DGG) / (GR3 * GR3);
        Real WDXG = -(Q1DX * GR3DG + Q1 * (GR3DXG - 2.0_rt * GR3DX * GR3DG / GR3)) / (GR3 * GR3);
        FSCR = -GAMI * Finf * W;
        Real FDX = -GAMI * (FinfDX * W + Finf * WDX);
        Real FDXX = -GAMI * (FinfDXX * W + 2.0_rt * FinfDX * WDX + Finf * WDXX);
        Real FDG = -Finf * W - GAMI * Finf * WDG;
        Real FDGG = -2.0_rt * Finf * WDG - GAMI * Finf * WDGG;
        if (std::abs(FDGG) < TINY) {
            FDGG = 0.0_rt; // 10.08.16: roundoff err.safeguard
        }
        Real FDXG = -FinfDX * W - Finf * WDX - GAMI * (FinfDX * WDG + Finf * WDXG);
        S_SCR = -GAMI * GAMI * Finf * WDG;
        USCR = S_SCR + FSCR;
        CVSCR = -GAMI * GAMI * FDGG;
        PSCR = (XSR * FDX + GAMI * FDG) / 3.0_rt;
        PDTSCR = GAMI * GAMI * (XSR * Finf * (FinfX * WDG + WDXG) - FDGG) / 3.0_rt;
        PDRSCR = (12.0_rt * PSCR + XSR * XSR * FDXX + 2.0_rt * XSR * GAMI * FDXG +
                  GAMI * GAMI * FDGG) / 9.0_rt;
    }

    void anharm8 (double GAMI, double TPT,
                  double& Fah, double& Uah, double& Pah,
                  double& CVah, double& PDTah, double& PDRah)
    {
        // ANHARMONIC free energy
        // Version 27.07.07
        // cleaned 16.06.09
        // Stems from ANHARM8b. Difference: AC = 0., B1 = .12 (.1217 - over accuracy)
        // Input: GAMI - ionic Gamma, TPT = Tp/T - ionic quantum parameter
        // Output: anharm.free en. Fah = F_{AH}/(N_i kT), internal energy Uah,
        //   pressure Pah = P_{AH}/(n_i kT), specific heat CVah  =  C_{V,AH}/(N_i k),
        //   PDTah  =  Pah  +  d Pah / d ln T, PDRah  =  Pah  +  d Pah / d ln\rho

        const int NM = 3;
        const Real AA[NM] = {10.9_rt, 247.0_rt, 1.765e5_rt}; // Farouki & Hamaguchi'93
        const Real B1 = 0.12_rt; // coeff.at \eta^2/\Gamma at T = 0

        Real CK = B1 / AA[0]; // fit coefficient
        Real TPT2 = TPT * TPT;
        Real TPT4 = TPT2 * TPT2;
        Real TQ = B1 * TPT2 / GAMI; // quantum dependence
        Real TK2 = CK * TPT2;
        Real SUP = std::exp(-TK2); // suppress.factor of class.anharmonicity

        Fah = 0.0_rt;
        Uah = 0.0_rt;
        Pah = 0.0_rt;
        CVah = 0.0_rt;
        PDTah = 0.0_rt;
        PDRah = 0.0_rt;

        Real SUPGN = SUP;
        for (int N = 1; N <= NM; ++N) {
            Real CN = (Real) N;
            SUPGN = SUPGN / GAMI; // SUP/Gamma^n
            Real ACN = AA[N-1];
            Fah = Fah - ACN / CN * SUPGN;
            Uah = Uah + (ACN * (1.0_rt + 2.0_rt * TK2 / CN)) * SUPGN;
            Real PN = AA[N-1] / 3.0_rt + TK2 * AA[N-1] / CN;
            Pah = Pah + PN * SUPGN;
            CVah = CVah + ((CN + 1.0_rt) * AA[N-1] + (4.0_rt - 2.0_rt / CN) * AA[N-1] * TK2 +
                           4.0_rt * AA[N-1] * CK * CK / CN * TPT4) * SUPGN;
            PDTah = PDTah + (PN * (1.0_rt + CN + 2.0_rt * TK2) - 2.0_rt / CN * AA[N-1] * TK2) * SUPGN;
            PDRah = PDRah + (PN * (1.0_rt - CN / 3.0_rt - TK2) + AA[N-1] / CN * TK2) * SUPGN;
        }

        Fah = Fah - TQ;
        Uah = Uah - TQ;
        Pah = Pah - TQ / 1.5_rt;
        PDRah = PDRah - TQ / 4.5_rt;
    }

    void hlfit12 (Real eta,
                  Real& F, Real& U, Real& CV, Real& S,
                  Real& U1, Real& CW, int LATTICE)
    {
        // Version 24.04.12
        // Stems from HLfit8 v.03.12.08;
        //   differences: E0 excluded from  U and F;
        //   U1 and d(CV)/d\ln(T) are added on the output.
        // Fit to thermal part of the thermodynamic functions.
        // Baiko, Potekhin, & Yakovlev (2001).
        // Zero-point lattice quantum energy 1.5u_1\eta EXCLUDED (unlike HLfit8).
        // Input: eta = Tp/T, LATTICE = 1 for bcc, 2 for fcc
        // Output: F and U (normalized to NkT) - due to phonon excitations,
        //   CV and S (normalized to Nk) in the HL model,
        //   U1 - the 1st phonon moment,
        //   CW = d(CV)/d\ln(T)

        const Real EPS = 1.e-5_rt;
        const Real TINY = 1.e-99_rt;

        Real CLM, ALPHA, BETA, GAMMA;
        Real A1, A2, A3, A4, A6, A8;
        Real B0, B2, B4, B5, B6, B7, C9, C11;

        if (LATTICE == 1) { // bcc lattice
            CLM = -2.49389_rt; // 3 * ln<\omega/\omega_p>
            U1 = 0.5113875_rt;
            ALPHA = 0.265764_rt;
            BETA = 0.334547_rt;
            GAMMA = 0.932446_rt;
            A1 = 0.1839_rt;
            A2 = 0.593586_rt;
            A3 = 0.0054814_rt;
            A4 = 5.01813e-4_rt;
            A6 = 3.9247e-7_rt;
            A8 = 5.8356e-11_rt;
            B0 = 261.66_rt;
            B2 = 7.07997_rt;
            B4 = 0.0409484_rt;
            B5 = 0.000397355_rt;
            B6 = 5.11148e-5_rt;
            B7 = 2.19749e-6_rt;
            C9 = 0.004757014_rt;
            C11 = 0.0047770935_rt;
        }
        else if (LATTICE == 2) { // fcc lattice
            CLM = -2.45373_rt;
            U1 = 0.513194_rt;
            ALPHA = 0.257591_rt;
            BETA = 0.365284_rt;
            GAMMA = 0.9167070_rt;
            A1 = 0.0_rt;
            A2 = 0.532535_rt;
            A3 = 0.0_rt;
            A4 = 3.76545e-4_rt;
            A6 = 2.63013e-7_rt;
            A8 = 6.6318e-11_rt;
            B0 = 303.20_rt;
            B2 = 7.7255_rt;
            B4 = 0.0439597_rt;
            B5 = 0.000114295_rt;
            B6 = 5.63434e-5_rt;
            B7 = 1.36488e-6_rt;
            C9 = 0.00492387_rt;
            C11 = 0.00437506_rt;
        }
        else {
            printf("HLfit: unknown lattice type\n");
            exit(1);
        }

        if (eta > 1.0_rt / EPS) { // asymptote of Eq.(13) of BPY'01
            U = 3.0_rt / (C11 * eta * eta * eta);
            F = -U / 3.0_rt;
            CV = 4.0_rt * U;
            S = U - F;
            return;
        }
        else if (eta < EPS) { // Eq.(17) of BPY'01
            if (eta < TINY) {
                printf("HLfit: eta is too small\n");
                exit(1);
            }
            F = 3.0_rt * std::log(eta) + CLM - 1.5_rt * U1 * eta + eta * eta / 24.0_rt;
            U = 3.0_rt - 1.5_rt * U1 * eta + eta * eta / 12.0_rt;
            CV = 3.0_rt - eta * eta / 12.0_rt;
            S = U - F;
            return;
        }

        Real eta2 = eta * eta;
        Real eta3 = eta2 * eta;
        Real eta4 = eta3 * eta;
        Real eta5 = eta4 * eta;
        Real eta6 = eta5 * eta;
        Real eta7 = eta6 * eta;
        Real eta8 = eta7 * eta;
        Real B9 = A6 * C9;
        Real B11 = A8 * C11;
        Real UP = 1.0_rt + A1 * eta + A2 * eta2 + A3 * eta3 + A4 * eta4 + A6 * eta6 + A8 * eta8;
        Real DN = B0 + B2 * eta2 + B4 * eta4 + B5 * eta5 + B6 * eta6 +
                  B7 * eta7 + eta8 * (B9 * eta + B11 * eta3);
        Real EA = std::exp(-ALPHA * eta);
        Real EB = std::exp(-BETA * eta);
        Real EG = std::exp(-GAMMA * eta);
        F = std::log(1.0_rt - EA) + std::log(1.0_rt - EB) + std::log(1.0_rt - EG) - UP / DN; // F_{thermal}/NT
        Real UP1 = A1 + 2.0_rt * A2 * eta + 3.0_rt * A3 * eta2 + 4.0_rt * A4 * eta3 +
                   6.0_rt * A6 * eta5 + 8. * A8 * eta7;
        Real UP2 = 2.0_rt * A2 + 6.0_rt * A3 * eta + 12.0_rt * A4 * eta2 + 30.0_rt * A6 * eta4 + 56.0_rt * A8 * eta6;
        Real UP3 = 6.0_rt * A3 + 24.0_rt * A4 * eta + 120.0_rt * A6 * eta3 + 336.0_rt * A8 * eta5;
        Real DN1 = 2.0_rt * B2 * eta + 4.0_rt * B4 * eta3 + 5.0_rt * B5 * eta4 + 6.0_rt * B6 * eta5 +
                   7.0_rt * B7 * eta6 + eta8 * (9.0_rt * B9 + 11.0_rt * B11 * eta2);
        Real DN2 = 2.0_rt * B2 + 12.0_rt * B4 * eta2 + 20. * B5 * eta3 + 30.0_rt * B6 * eta4 +
                   42.0_rt * B7 * eta5 + 72.0_rt * B9 * eta7 + 110.0_rt * B11 * eta8 * eta;
        Real DN3 = 24.0_rt * B4 * eta + 60.0_rt * B5 * eta2 + 120.0_rt * B6 * eta3 +
                   210.0_rt * B7 * eta4 + 504.0_rt * B9 * eta6 + 990.0_rt * B11 * eta8;
        Real DF1 = ALPHA * EA / (1.0_rt - EA) + BETA * EB / (1.0_rt - EB) + GAMMA * EG / (1.0_rt - EG) -
                   (UP1 * DN - DN1 * UP) / (DN * DN); // int.en./NT/eta  =  df/d\eta
        Real DF2 = ALPHA * ALPHA * EA / ((1.0_rt - EA) * (1.0_rt - EA)) + BETA * BETA * EB /
                   ((1.0_rt - EB) * (1.0_rt - EB)) + GAMMA * GAMMA * EG / ((1.0_rt - EG) * (1.0_rt - EG)) +
                   ((UP2 * DN - DN2 * UP) * DN - 2.0_rt * (UP1 * DN - DN1 * UP) * DN1) / (DN * DN * DN); // -d2f/d\eta^2
        U = DF1 * eta;
        CV = DF2 * eta2;
        Real DF3 = -ALPHA * ALPHA * ALPHA * EA / std::pow(1.0_rt - EA, 3) * (1.0_rt + EA) -
                   BETA * BETA * BETA * EB / std::pow(1.0_rt - EB, 3) * (1.0_rt + EB) -
                   GAMMA * GAMMA * GAMMA * EG / std::pow(1.0_rt - EG, 3) * (1.0_rt + EG) +
                   UP3 / DN - (3.0_rt * UP2 * DN1 + 3.0_rt * UP1 * DN2 + UP * DN3) / (DN * DN) +
                   6.0_rt * DN1 * (UP1 * DN1 + UP * DN2) / (DN * DN * DN) -
                   6.0_rt * UP * DN1 * DN1 * DN1 / (DN * DN * DN * DN); // -d3f/d\eta^3
        CW = -2.0_rt * CV - eta3 * DF3;
        S = U - F;
    }

    /*
      subroutine FHARM12(GAMI,TPT, &
        Fharm,Uharm,Pharm,CVth,Sth,PDTharm,PDRharm)
// Thermodynamic functions of a harmonic crystal, incl.stat.Coul.lattice
//
//                                                       Version 27.04.12
// Stems from FHARM8 v.15.02.08
// Replaced HLfit8 with HLfit12: rearranged output.
// Input: GAMI - ionic Gamma, TPT = T_{p,i}/T
// Output: Fharm = F/(N_i T), Uharm = U/(N_i T), Pharm = P/(n_i T),
// CVth = C_V/N_i, Sharm = S/N_i
// PDTharm  =  Pharm  +  d Pharm / d ln T, PDRharm  =  Pharm  +  d Pharm/d ln\rho
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(CM = .895929256d0) // Madelung
      call HLfit12(TPT,F,U,CVth,Sth,U1,CW,1)
      U0 = -CM * GAMI // perfect lattice
      E0 = 1.5d0 * U1 * TPT // zero-point energy
      Uth = U + E0
      Fth = F + E0
      Uharm = U0 + Uth
      Fharm = U0 + Fth
      Pharm = U0 / 3.0_rt + Uth / 2.0_rt
      PDTharm = 0.5_rt * CVth
      PDRharm = U0 / 2.25d0 + .75d0 * Uth - .25d0 * CVth
      return
      end
    */

    void cormix (Real RS, Real GAME, Real Zmean, Real Z2mean, Real Z52, Real Z53, Real Z321,
                 Real& FMIX, Real& UMIX, Real& PMIX, Real& CVMIX, Real& PDTMIX, Real& PDRMIX)
    {
        // Version 02.07.09
        // Correction to the linear mixing rule for moderate to small Gamma
        // Input: RS = r_s (if RS = 0, then OCP, otherwise EIP)
        //        GAME = \Gamma_e
        //        Zmean = <Z> (average Z of all ions, without electrons)
        //        Z2mean = <Z^2>, Z52 = <Z^2.5>, Z53 = <Z^{5/3}>, Z321 = <Z(Z + 1)^1.5>
        // Output: FMIX = \Delta f - corr.to the reduced free energy f = F/N_{ion}kT
        //         UMIX = \Delta u - corr.to the reduced internal energy u
        //         PMIX = \Delta u - corr.to the reduced pressure P = P/n_{ion}kT
        //         CVMIX = \Delta c - corr.to the reduced heat capacity c_V
        //         PDTMIX = (1/n_{ion}kT)d\Delta P / d ln T
        //                =  \Delta p  +   d \Delta p / d ln T
        //         PDRMIX = (1/n_{ion}kT)d\Delta P / d ln n_e
        // (composition is assumed fixed: Zmean,Z2mean,Z52,Z53 = constant)

        const Real TINY = 1.e-9_rt;
        Real GAMImean = GAME * Z53;

        Real Dif0;
        if (RS < TINY) { // OCP
            Dif0 = Z52 - std::sqrt(Z2mean * Z2mean * Z2mean / Zmean);
        }
        else {
            Dif0 = Z321 - std::sqrt(std::pow(Z2mean + Zmean, 3) / Zmean);
        }

        Real DifR = Dif0 / Z52;
        Real DifFDH = Dif0 * GAME * std::sqrt(GAME / 3.0_rt); // F_DH - F_LM(DH)
        Real D = Z2mean / (Zmean * Zmean);
        if (std::abs(D - 1.0_rt) < TINY) { // no correction
            FMIX = 0.0_rt;
            UMIX = 0.0_rt;
            PMIX = 0.0_rt;
            CVMIX = 0.0_rt;
            PDTMIX = 0.0_rt;
            PDRMIX = 0.0_rt;
            return;
        }

        Real P3 = std::pow(D, -0.2_rt);
        Real D0 = (2.6_rt * DifR + 14.0_rt * DifR * DifR * DifR) / (1.0_rt - P3);
        Real GP = D0 * std::pow(GAMImean, P3);
        Real FMIX0 = DifFDH / (1.0_rt + GP);
        Real Q = D * D * 0.0117_rt;
        Real R = 1.5_rt / P3 - 1.0_rt;
        Real GQ = Q * GP;
        FMIX = FMIX0 / std::pow(1.0_rt + GQ, R);
        Real G = 1.5_rt - P3 * GP / (1.0_rt + GP) - R * P3 * GQ / (1.0_rt + GQ);
        UMIX = FMIX * G;
        PMIX = UMIX / 3.0_rt;
        Real GDG = -P3 * P3 * (GP / ((1.0_rt + GP) * (1.0_rt + GP)) + R * GQ / ((1.0_rt + GQ) * (1.0_rt + GQ))); // d G /d ln Gamma
        Real UDG = UMIX * G + FMIX * GDG; // d u_mix /d ln Gamma
        CVMIX = UMIX - UDG;
        PDTMIX = PMIX - UDG / 3.0_rt;
        PDRMIX = PMIX + UDG / 9.0_rt;
    }
}
