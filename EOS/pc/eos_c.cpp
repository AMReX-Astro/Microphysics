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
}
