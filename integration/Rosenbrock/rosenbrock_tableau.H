#ifndef ROSENBROCK_TABLEAU_H
#define ROSENBROCK_TABLEAU_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace rosenbrock {

struct ros2s_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.292893218813452_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 2 ? 0.585786437626905_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i == 2 && j == 1) {
            return 2.0000000000000036_rt;
        }
        if (i == 3 && j == 1) {
            return 6.828427124746214_rt;
        }
        if (i == 3 && j == 2) {
            return 3.4142135623731007_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -6.828427124746214_rt;
        }
        if (i == 3 && j == 1) {
            return -10.949747468305889_rt;
        }
        if (i == 3 && j == 2) {
            return -7.535533905932761_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 6.828427124746214_rt;
        }
        if (i == 2) {
            return 3.414213562373101_rt;
        }
        return 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        if (i == 1) {
            return -0.23570226039551292_rt;
        }
        if (i == 2) {
            return -0.23570226039551567_rt;
        }
        return -0.13807118745769906_rt;
    }
};

struct ros2_tableau {
    static constexpr int stages = 2;
    static constexpr amrex::Real gamma = 1.7071067811865475_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        (void) i;
        return 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i == 2 && j == 1) {
            return 0.585786437626905_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -1.1715728752538102_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 0.8786796564403575_rt;
        }
        return 0.2928932188134525_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        (void) i;
        return 0.2928932188134525_rt;
    }
};

struct sandu_a_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.7886751345948129_rt; // (3 + sqrt(3)) / 6

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 1 ? 0.0_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if ((i == 2 || i == 3) && j == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return 0.9282032302755092_rt; // -6 + 4 * sqrt(3)
        }
        if (i == 3 && (j == 1 || j == 2)) {
            return -0.4641016151377544_rt; // 3 - 2 * sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        if (i == 3) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        return i == 3 ? 1.0_rt : 0.0_rt;
    }
};

struct sandu_b_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.7886751345948129_rt; // (3 + sqrt(3)) / 6

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 1 ? 0.0_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if ((i == 2 || i == 3) && j == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 3 && j == 1) {
            return -0.5358983848622456_rt; // -4 + 2 * sqrt(3)
        }
        if (i == 3 && j == 2) {
            return -0.7320508075688772_rt; // 1 - sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 1.2679491924311228_rt; // 3 - sqrt(3)
        }
        if (i == 3) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        return i == 3 ? 1.0_rt : 0.0_rt;
    }
};

struct sandu_c_tableau {
    static constexpr int stages = 3;
    static constexpr amrex::Real gamma = 0.78867513459481275_rt; // (3 + sqrt(3)) / 6

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        if (i == 3) {
            return 0.6666666666666666_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i == 3 && j == 1) {
            return -0.4335531486700170_rt; // (1272 - 823 * sqrt(3)) / 354
        }
        if (i == 3 && j == 2) {
            return 1.7222282832648295_rt; // 3 * (-51 + 49 * sqrt(3)) / 59
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -0.3264246859454366_rt; // -1 + 7 * sqrt(3) / 18
        }
        if (i == 3 && j == 1) {
            return 0.4138899169340995_rt; // (25 - 13 * sqrt(3)) / 6
        }
        if (i == 3 && j == 2) {
            return 1.6076951545867360_rt; // 12 - 6 * sqrt(3)
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return -7.1285171234651825_rt; // (-12089 + 5037 * sqrt(3)) / 472
        }
        if (i == 2) {
            return 8.4030361763035106_rt; // 9 * (344 - 135 * sqrt(3)) / 118
        }
        return 0.9509618943233421_rt; // 3 * (3 - sqrt(3)) / 4
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        if (i == 1) {
            return -8.396466315896305_rt;
        }
        if (i == 2) {
            return 8.4030361763035106_rt;
        }
        return 0.9509618943233421_rt;
    }
};

struct sandu_d_tableau {
    static constexpr int stages = 4;
    static constexpr amrex::Real gamma = 0.5_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        return i == 1 ? 0.0_rt : 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i > 1 && j == 1) {
            return 2.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -1.3333333333333333_rt; // -4 / 3
        }
        if (i == 3 && j == 1) {
            return -3.3333333333333333_rt; // -10 / 3
        }
        if (i == 3 && j == 2) {
            return -2.0_rt;
        }
        if (i == 4 && j == 1) {
            return -0.5_rt;
        }
        if (i == 4 && j == 3) {
            return 1.5_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return 2.0_rt;
        }
        if (i == 4) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        if (i == 1) {
            return -0.6666666666666666_rt; // 2 - 8 / 3
        }
        if (i == 2 || i == 3) {
            return -1.0_rt;
        }
        return 1.3333333333333333_rt; // 1 - (-1 / 3)
    }
};

struct rodas5p_tableau {
    static constexpr int stages = 8;
    static constexpr amrex::Real gamma = 0.21193756319429014_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        if (i == 2) {
            return 0.6358126895828704_rt;
        }
        if (i == 3) {
            return 0.4095798393397535_rt;
        }
        if (i == 4) {
            return 0.9769306725060716_rt;
        }
        if (i == 5) {
            return 0.4288403609558664_rt;
        }
        if (i >= 6) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        if (i == 2 && j == 1) {
            return 3.0_rt;
        }
        if (i == 3 && j == 1) {
            return 2.849394379747939_rt;
        }
        if (i == 3 && j == 2) {
            return 0.45842242204463923_rt;
        }
        if (i == 4 && j == 1) {
            return -6.954028509809101_rt;
        }
        if (i == 4 && j == 2) {
            return 2.489845061869568_rt;
        }
        if (i == 4 && j == 3) {
            return -10.358996098473584_rt;
        }
        if (i == 5 && j == 1) {
            return 2.8029986275628964_rt;
        }
        if (i == 5 && j == 2) {
            return 0.5072464736228206_rt;
        }
        if (i == 5 && j == 3) {
            return -0.3988312541770524_rt;
        }
        if (i == 5 && j == 4) {
            return -0.04721187230404641_rt;
        }
        if ((i == 6 || i == 7 || i == 8) && j == 1) {
            return -7.502846399306121_rt;
        }
        if ((i == 6 || i == 7 || i == 8) && j == 2) {
            return 2.561846144803919_rt;
        }
        if ((i == 6 || i == 7 || i == 8) && j == 3) {
            return -11.627539656261098_rt;
        }
        if ((i == 6 || i == 7 || i == 8) && j == 4) {
            return -0.18268767659942256_rt;
        }
        if ((i == 6 || i == 7 || i == 8) && j == 5) {
            return 0.030198172008377946_rt;
        }
        if ((i == 7 || i == 8) && j == 6) {
            return 1.0_rt;
        }
        if (i == 8 && j == 7) {
            return 1.0_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        if (i == 2 && j == 1) {
            return -14.155112264123755_rt;
        }
        if (i == 3 && j == 1) {
            return -17.97296035885952_rt;
        }
        if (i == 3 && j == 2) {
            return -2.859693295451294_rt;
        }
        if (i == 4 && j == 1) {
            return 147.12150275711716_rt;
        }
        if (i == 4 && j == 2) {
            return -1.41221402718213_rt;
        }
        if (i == 4 && j == 3) {
            return 71.68940251302358_rt;
        }
        if (i == 5 && j == 1) {
            return 165.43517024871676_rt;
        }
        if (i == 5 && j == 2) {
            return -0.4592823456491126_rt;
        }
        if (i == 5 && j == 3) {
            return 42.90938336958603_rt;
        }
        if (i == 5 && j == 4) {
            return -5.961986721573306_rt;
        }
        if (i == 6 && j == 1) {
            return 24.854864614690072_rt;
        }
        if (i == 6 && j == 2) {
            return -3.0009227002832186_rt;
        }
        if (i == 6 && j == 3) {
            return 47.4931110020768_rt;
        }
        if (i == 6 && j == 4) {
            return 5.5814197821558125_rt;
        }
        if (i == 6 && j == 5) {
            return -0.6610691825249471_rt;
        }
        if (i == 7 && j == 1) {
            return 30.91273214028599_rt;
        }
        if (i == 7 && j == 2) {
            return -3.1208243349937974_rt;
        }
        if (i == 7 && j == 3) {
            return 77.79954646070892_rt;
        }
        if (i == 7 && j == 4) {
            return 34.28646028294783_rt;
        }
        if (i == 7 && j == 5) {
            return -19.097331116725623_rt;
        }
        if (i == 7 && j == 6) {
            return -28.087943162872662_rt;
        }
        if (i == 8 && j == 1) {
            return 37.80277123390563_rt;
        }
        if (i == 8 && j == 2) {
            return -3.2571969029072276_rt;
        }
        if (i == 8 && j == 3) {
            return 112.26918849496327_rt;
        }
        if (i == 8 && j == 4) {
            return 66.9347231244047_rt;
        }
        if (i == 8 && j == 5) {
            return -40.06618937091002_rt;
        }
        if (i == 8 && j == 6) {
            return -54.66780262877968_rt;
        }
        if (i == 8 && j == 7) {
            return -9.48861652309627_rt;
        }
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        if (i == 1) {
            return -7.502846399306121_rt;
        }
        if (i == 2) {
            return 2.561846144803919_rt;
        }
        if (i == 3) {
            return -11.627539656261098_rt;
        }
        if (i == 4) {
            return -0.18268767659942256_rt;
        }
        if (i == 5) {
            return 0.030198172008377946_rt;
        }
        return 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        return i == 8 ? 1.0_rt : 0.0_rt;
    }
};

struct rosenbrock_euler_tableau {
    static constexpr int stages = 1;
    static constexpr amrex::Real gamma = 1.0_rt;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real ct (const int i) {
        (void) i;
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real a (const int i, const int j) {
        (void) i;
        (void) j;
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real c (const int i, const int j) {
        (void) i;
        (void) j;
        return 0.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real b (const int i) {
        (void) i;
        return 1.0_rt;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static constexpr amrex::Real e (const int i) {
        (void) i;
        return 1.0_rt;
    }
};

using coefficients = ros2s_tableau;

}

#endif
