#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

extern "C"
{
    // Inverse Fermi integrals
    void ferinv7 (double F, int N,
                  double& X,
                  double& XDF,
                  double& XDFF)
    {
        // Version 24.05.07
        // X_q(f)=F^{-1}_q(f) : H.M.Antia 93 ApJS 84, 101
        // q=N-1/2=-1/2,1/2,3/2,5/2 (N=0,1,2,3)
        // Input: F - argument, N=q+1/2
        // Output: X=X_q, XDF=dX/df, XDFF=d^2 X / df^2
        // Relative error:  N = 0     1       2       3
        //        for X:    3.e-9,  4.2e-9, 2.3e-9, 6.2e-9
        // jump at f=4:
        //        for XDF:  6.e-7,  5.4e-7, 9.6e-8, 3.1e-7
        //        for XDFF: 4.7e-5, 4.8e-5, 2.3e-6, 1.5e-6

        const double A[4][6] = {{-1.570044577033e4,  1.001958278442e4, -2.805343454951e3,
                                  4.121170498099e2, -3.174780572961e1,  1.e0}, // X_{-1/2}
                                { 1.999266880833e4,  5.702479099336e3,  6.610132843877e2,
                                  3.818838129486e1,  1.e0,              0.0}, // X_{1/2}
                                { 1.715627994191e2,  1.125926232897e2,  2.056296753055e1,
                                  1.e0,              0.0,               0.0},
                                { 2.138969250409e2,  3.539903493971e1,  1.e0,
                                  0.0,               0.0,               0.0}}; // X_{5/2}

        const double B[4][7] = {{-2.782831558471e4,  2.886114034012e4, -1.274243093149e4,
                                  3.063252215963e3, -4.225615045074e2,  3.168918168284e1,
                                 -1.008561571363e0}, // X_{-1/2}
                                { 1.771804140488e4, -2.014785161019e3,  9.130355392717e1,
                                 -1.670718177489e0,  0.0,               0.0,
                                  0.0}, // X_{1/2}
                                { 2.280653583157e2,  1.193456203021e2,  1.16774311354e1,
                                 -3.226808804038e-1, 3.519268762788e-3, 0.0,
                                  0.0}, // X_{3/2}
                                { 7.10854551271e2,   9.873746988121e1,  1.067755522895e0,
                                 -1.182798726503e-2, 0.0,               0.0,
                                  0.0}}; // X_{5/2}

        const double C[4][7] = {{ 2.206779160034e-8, -1.437701234283e-6,  6.103116850636e-5,
                                 -1.169411057416e-3,  1.814141021608e-2, -9.588603457639e-2,
                                  1.e0},
                                {-1.277060388085e-2,  7.187946804945e-2, -4.262314235106e-1,
                                  4.997559426872e-1, -1.285579118012e0,  -3.930805454272e-1,
                                  1.e0},
                                {-6.321828169799e-3, -2.183147266896e-2, -1.05756279932e-1,
                                 -4.657944387545e-1, -5.951932864088e-1,  3.6844711771e-1,
                                  1.e0},
                                {-3.312041011227e-2,  1.315763372315e-1, -4.820942898296e-1,
                                  5.099038074944e-1,  5.49561349863e-1,  -1.498867562255e0,
                                  1.e0}};

        const double D[4][7] = {{ 8.827116613576e-8, -5.750804196059e-6,  2.429627688357e-4,
                                 -4.601959491394e-3,  6.932122275919e-2, -3.217372489776e-1,
                                  3.124344749296e0}, // X_{-1/2}
                                {-9.745794806288e-3,  5.485432756838e-2, -3.29946624326e-1,
                                  4.077841975923e-1, -1.145531476975e0,  -6.067091689181e-2,
                                  0.0},
                                {-4.381942605018e-3, -1.5132365041e-2,   -7.850001283886e-2,
                                 -3.407561772612e-1, -5.074812565486e-1, -1.387107009074e-1,
                                  0.0},
                                {-2.315515517515e-2,  9.198776585252e-2, -3.835879295548e-1,
                                  5.415026856351e-1, -3.847241692193e-1,  3.739781456585e-2,
                                 -3.008504449098e-2}}; // X_{5/2}

        const int LA[4] = {5, 4, 3, 2};
        const int LB[4] = {6, 3, 4, 3};
        const int LD[4] = {6, 5, 5, 6};

         if (N < 0 || N > 3) {
             printf("FERINV7: Invalid subscript\n");
             exit(1);
         }
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
             for (int i = LA[N]; i >= 0; --i) {
                 UP = UP * T + A[N][i];
                 if (i >= 1) {
                     UP1 = UP1 * T + A[N][i] * i;
                 }
                 if (i >= 2) {
                     UP2 = UP2 * T + A[N][i] * i * (i-1);
                 }
             }
             for (int i = LB[N]; i >= 0; --i) {
                 DOWN = DOWN * T + B[N][i];
                 if (i >= 1) {
                     DOWN1 = DOWN1 * T + B[N][i] * i;
                 }
                 if (i >= 2) {
                     DOWN2 = DOWN2 * T + B[N][i] * i * (i-1);
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
                 UP = UP * T + C[N][i];
                 if (i >= 1) {
                     UP1 = UP1 * T + C[N][i] * i;
                 }
                 if (i >= 2) {
                     UP2 = UP2 * T + C[N][i] * i * (i-1);
                 }
             }
             for (int i = LD[N]; i >= 0; --i) {
                 DOWN = DOWN * T + D[N][i];
                 if (i >= 1) {
                     DOWN1 = DOWN1 * T + D[N][i] * i;
                 }
                 if (i >= 2) {
                     DOWN2 = DOWN2 * T + D[N][i] * i * (i-1);
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
