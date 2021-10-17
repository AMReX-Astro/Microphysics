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

    void blin9a (double TEMP, double CHI,
                 double& W0, double& W0DX, double& W0DT, double& W0DXX,
                 double& W0DTT, double& W0DXT,
                 double& W1, double& W1DX, double& W1DT, double& W1DXX,
                 double& W1DTT, double& W1DXT,
                 double& W2, double& W2DX, double& W2DT, double& W2DXX,
                 double& W2DTT, double& W2DXT,
                 double& W0XXX, double& W0XTT, double& W0XXT)
    {
        // Version 19.01.10
        // First part of blin9: small CHI. Stems from blin9 v.24.12.08
        const double AC[3][5] = {{0.37045057_rt,  0.41258437_rt,
                                  9.777982e-2_rt, 5.3734153e-3_rt, 3.8746281e-5_rt},  // c_i^0
                                 {0.39603109_rt,  0.69468795_rt,
                                  0.22322760_rt,  1.5262934e-2_rt, 1.3081939e-4_rt},  // c_i^1
                                 {0.76934619_rt,  1.7891437_rt,
                                  0.70754974_rt,  5.6755672e-2_rt, 5.5571480e-4_rt}}; // c_i^2

        const double AU[3][5] = {{0.43139881_rt,  1.7597537_rt,
                                  4.10446540_rt,  7.7467038_rt, 13.457678_rt},  // \chi_i^0
                                 {0.81763176_rt,  2.4723339_rt,
                                  5.11600610_rt,  9.0441465_rt, 15.049882_rt},  // \chi_i^1
                                 {1.25584610_rt,  3.2070406_rt,
                                  6.12390820_rt, 10.3161260_rt, 16.597079_rt}}; // \chi_i^2

        const double AA[3][5] = {{std::exp(-AU[0][0]), std::exp(-AU[0][1]),
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

    void fermi10 (double X, double XMAX, double& FP, double& FM)
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
}
