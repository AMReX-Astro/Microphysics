#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

extern "C"
{
    // Inverse Fermi integrals with q=1/2
    void ferinv7 (double F,
                  double& X,
                  double& XDF,
                  double& XDFF)
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

        const double A[6] = { 1.999266880833e4,  5.702479099336e3,  6.610132843877e2,
                              3.818838129486e1,  1.e0,              0.0};
        const double B[7] = { 1.771804140488e4, -2.014785161019e3,  9.130355392717e1,
                             -1.670718177489e0,  0.0,               0.0,
                              0.0};
        const double C[7] = {-1.277060388085e-2,  7.187946804945e-2, -4.262314235106e-1,
                              4.997559426872e-1, -1.285579118012e0,  -3.930805454272e-1,
                              1.e0};
        const double D[7] = {-9.745794806288e-3,  5.485432756838e-2, -3.29946624326e-1,
                              4.077841975923e-1, -1.145531476975e0,  -6.067091689181e-2,
                              0.0};
        const int LA = 4;
        const int LB = 3;
        const int LD = 5;

        const int N = 1;

        if (F <= 0.0) {
            printf("FERINV7: Non-positive argument\n");
            exit(1);
        }
        if (F < 4.0) {
            double T = F;
            double UP = 0.0;
            double UP1 = 0.0;
            double UP2 = 0.0;
            double DOWN = 0.0;
            double DOWN1 = 0.0;
            double DOWN2 = 0.0;
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
            XDF = 1.0 / T + UP1 / UP - DOWN1 / DOWN;
            XDFF = -1.0 / (T * T) + UP2 / UP - (UP1 / UP) * (UP1 / UP) -
                DOWN2 / DOWN + (DOWN1 / DOWN) * (DOWN1 / DOWN);
        }
        else {
            double P = -1.0 / (0.5 + (double) N); // = -1/(1+\nu) = power index
            double T = std::pow(F, P); // t  - argument of the rational fraction
            double T1 = P * T / F; // dt/df
            double T2 = P * (P - 1.0) * T / (F * F); // d^2 t / df^2
            double UP = 0.0;
            double UP1 = 0.0;
            double UP2 = 0.0;
            double DOWN = 0.0;
            double DOWN1 = 0.0;
            double DOWN2 = 0.0;
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
            double R = UP / DOWN;
            double R1 = (UP1 - UP * DOWN1 / DOWN) / DOWN; // dR/dt
            double R2 = (UP2 - (2.0 * UP1 * DOWN1 + UP * DOWN2) / DOWN +
                         2.0 * UP * (DOWN1 / DOWN) * (DOWN1 / DOWN)) / DOWN;
            X = R/T;
            double RT = (R1 - R / T) / T;
            XDF = T1 * RT;
            XDFF = T2 * RT + T1 * T1 * (R2 - 2.0 * RT) / T;
        }
    }

    void chemfit (double DENS, double TEMP, double& CHI)
    {
        // Version 29.08.15
        // Fit to the chemical potential of free electron gas described in:
        //     G.Chabrier & A.Y.Potekhin, Phys.Rev.E 58, 4941 (1998)
        // Stems from CHEMFIT v.10.10.96. The main difference - derivatives.
        // Input:  DENS - electron density [a.u.=6.7483346e24 cm^{-3}],
        // TEMP - temperature [a.u.=2Ryd=3.1577e5 K]
        // Output: CHI = CMU1 / TEMR, where CMU1 = \mu-1 - chem.pot.w/o rest-energy

        const double C13 = 1.0 / 3.0;
        const double PARA = 1.612;
        const double PARB = 6.192;
        const double PARC = 0.0944;
        const double PARF = 5.535;
        const double PARG = 0.698;
        const double XEPST = 228.0; // the largest argument of e^{-X}

        double DENR = DENS / 2.5733806e6; // n_e in rel.un.=\lambda_{Compton}^{-3}
        double TEMR = TEMP / 1.8778865e4; // T in rel.un.=(mc^2/k)=5.93e9 K

        double PF0 = std::pow(29.6088132 * DENR, C13); // Classical Fermi momentum
        double TF;
        if (PF0 > 1.e-4) {
            TF = std::sqrt(1.0 + PF0 * PF0) - 1.0; // Fermi temperature
        }
        else {
            TF = 0.50 * PF0 * PF0;
        }

        double THETA = TEMR / TF;
        double THETA32 = THETA * std::sqrt(THETA);
        double Q2 = 12.0 + 8.0 / THETA32;
        double T1 = 0.0;
        if (THETA < XEPST) {
            T1 = std::exp(-THETA);
        }
        double U3 = T1 * T1 + PARA;
        double THETAC = std::pow(THETA, PARC);
        double THETAG = std::pow(THETA, PARG);
        double D3 = PARB * THETAC * T1 * T1 + PARF * THETAG;
        double Q3 = 1.365568127 - U3 / D3; // 1.365...=2/\pi^{1/3}
        double Q1;
        if (THETA > 1.e-5) {
            Q1 = 1.5 * T1 / (1.0 - T1);
        }
        else {
            Q1 = 1.5 / THETA;
        }
        double SQT = std::sqrt(TEMR);
        double G = (1.0 + Q2 * TEMR * Q3 + Q1 * SQT) * TEMR;
        double H = (1.0 + 0.5 * TEMR / THETA) * (1.0 + Q2 * TEMR);
        double CT = 1.0 + G / H;
        double F = 2.0 * C13 / THETA32;
        double X, XDF, XDFF;
        ferinv7(F, X, XDF, XDFF);
        CHI = X                    // Non-relativistic result
            - 1.50 * std::log(CT); // Relativistic fit
    }
}
