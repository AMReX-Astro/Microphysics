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
    // Inverse Fermi integrals with q=1/2
    void ferinv7 (Real F,
                  Real& X,
                  Real& XDF,
                  Real& XDFF)
    {
        // Version 24.05.07
        // X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84, 101
        // Input: F - argument
        // Output: X=X_q, XDF=dX/df, XDFF=d^2 X / df^2
        // Relative error:
        //        for X:    4.2e-9
        // jump at f=4:
        //        for XDF:  5.4e-7
        //        for XDFF: 4.8e-5

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
            printf("FERINV7: Non-positive argument\n");
            exit(1);
        }
        if (F < 4.0_rt) {
            Real T = F;
            Real UP = 0.0_rt;
            Real UP1 = 0.0_rt;
            Real UP2 = 0.0_rt;
            Real DOWN = 0.0_rt;
            Real DOWN1 = 0.0_rt;
            Real DOWN2 = 0.0_rt;
            for (int i = LA; i >= 0; --i) {
                UP = UP * T + A[i];
                if (i >= 1) {
                    UP1 = UP1 * T + A[i] * i;
                }
                if (i >= 2) {
                    UP2 = UP2 * T + A[i] * i * (i-1);
                }
            }
            for (int i = LB; i >= 0; --i) {
                DOWN = DOWN * T + B[i];
                if (i >= 1) {
                    DOWN1 = DOWN1 * T + B[i] * i;
                }
                if (i >= 2) {
                    DOWN2 = DOWN2 * T + B[i] * i * (i-1);
                }
            }
            X = std::log(T * UP / DOWN);
            XDF = 1.0_rt / T + UP1 / UP - DOWN1 / DOWN;
            XDFF = -1.0_rt / (T * T) + UP2 / UP - (UP1 / UP) * (UP1 / UP) -
                DOWN2 / DOWN + (DOWN1 / DOWN) * (DOWN1 / DOWN);
        }
        else {
            Real P = -1.0_rt / (0.5_rt + (Real) N); // = -1/(1+\nu) = power index
            Real T = std::pow(F, P); // t  - argument of the rational fraction
            Real T1 = P * T / F; // dt/df
            Real T2 = P * (P - 1.0_rt) * T / (F * F); // d^2 t / df^2
            Real UP = 0.0_rt;
            Real UP1 = 0.0_rt;
            Real UP2 = 0.0_rt;
            Real DOWN = 0.0_rt;
            Real DOWN1 = 0.0_rt;
            Real DOWN2 = 0.0_rt;
            for (int i = 6; i >= 0; --i) {
                UP = UP * T + C[i];
                if (i >= 1) {
                    UP1 = UP1 * T + C[i] * i;
                }
                if (i >= 2) {
                    UP2 = UP2 * T + C[i] * i * (i-1);
                }
            }
            for (int i = LD; i >= 0; --i) {
                DOWN = DOWN * T + D[i];
                if (i >= 1) {
                    DOWN1 = DOWN1 * T + D[i] * i;
                }
                if (i >= 2) {
                    DOWN2 = DOWN2 * T + D[i] * i * (i-1);
                }
            }
            Real R = UP / DOWN;
            Real R1 = (UP1 - UP * DOWN1 / DOWN) / DOWN; // dR/dt
            Real R2 = (UP2 - (2.0_rt * UP1 * DOWN1 + UP * DOWN2) / DOWN +
                       2.0_rt * UP * (DOWN1 / DOWN) * (DOWN1 / DOWN)) / DOWN;
            X = R/T;
            Real RT = (R1 - R / T) / T;
            XDF = T1 * RT;
            XDFF = T2 * RT + T1 * T1 * (R2 - 2.0_rt * RT) / T;
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
        Real X, XDF, XDFF;
        ferinv7(F, X, XDF, XDFF);
        CHI = X                      // Non-relativistic result
            - 1.5_rt * std::log(CT); // Relativistic fit
    }
}
