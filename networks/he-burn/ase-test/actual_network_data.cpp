#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, 0, 2, -1, -1, 3, 59,  // p_C12_to_N13_reaclib
        -1, 1, 2, -1, -1, 4, 60,  // He4_C12_to_O16_reaclib
        -1, 1, 4, -1, -1, 6, 61,  // He4_O16_to_Ne20_reaclib
        -1, 1, 6, -1, -1, 8, 63,  // He4_Ne20_to_Mg24_reaclib
        -1, 0, 7, -1, -1, 8, 62,  // p_Na23_to_Mg24_reaclib
        -1, 1, 8, -1, -1, 10, 65,  // He4_Mg24_to_Si28_reaclib
        -1, 0, 9, -1, -1, 10, 64,  // p_Al27_to_Si28_reaclib
        -1, 1, 10, -1, -1, 12, 67,  // He4_Si28_to_S32_reaclib
        -1, 0, 11, -1, -1, 12, 66,  // p_P31_to_S32_reaclib
        -1, 2, 2, -1, 0, 7, 73,  // C12_C12_to_p_Na23_reaclib
        -1, 2, 2, -1, 1, 6, 72,  // C12_C12_to_He4_Ne20_reaclib
        -1, 1, 3, -1, 0, 4, 69,  // He4_N13_to_p_O16_reaclib
        -1, 2, 4, -1, 0, 9, 76,  // C12_O16_to_p_Al27_reaclib
        -1, 2, 4, -1, 1, 8, 75,  // C12_O16_to_He4_Mg24_reaclib
        -1, 4, 4, -1, 0, 11, 79,  // O16_O16_to_p_P31_reaclib
        -1, 4, 4, -1, 1, 10, 78,  // O16_O16_to_He4_Si28_reaclib
        -1, 0, 6, -1, 1, 5, -1,  // p_Ne20_to_He4_F17_reaclib
        -1, 0, 7, -1, 1, 6, 71,  // p_Na23_to_He4_Ne20_reaclib
        -1, 0, 9, -1, 1, 8, 74,  // p_Al27_to_He4_Mg24_reaclib
        -1, 0, 11, -1, 1, 10, 77,  // p_P31_to_He4_Si28_reaclib
        1, 1, 1, -1, -1, 2, 68,  // He4_He4_He4_to_C12_reaclib
        -1, -1, -1, -1, -1, -1, -1,  // C12_C12_to_n_Mg23_removed
        -1, -1, -1, -1, -1, -1, -1,  // O16_O16_to_n_S31_removed
        -1, -1, -1, -1, -1, -1, -1,  // C12_O16_to_n_Si27_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_S32_to_Ar36_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Cl35_to_Ar36_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Cl35_to_He4_S32_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Ar36_to_Ca40_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_K39_to_Ca40_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_K39_to_He4_Ar36_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Ca40_to_Ti44_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Sc43_to_Ti44_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Sc43_to_He4_Ca40_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Ti44_to_Cr48_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Ti44_to_p_V47_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_V47_to_Cr48_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Cr48_to_Fe52_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Cr48_to_p_Mn51_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Mn51_to_Fe52_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Fe52_to_Ni56_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Fe52_to_p_Co55_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Co55_to_Ni56_removed
        -1, 0, 4, -1, -1, 5, 83,  // O16_p_to_F17_iliadis
        -1, 1, 12, -1, -1, 13, 45,  // S32_He4_to_Ar36_approx
        -1, -1, 13, -1, 1, 12, -1,  // Ar36_to_S32_He4_approx
        -1, 1, 13, -1, -1, 14, 47,  // Ar36_He4_to_Ca40_approx
        -1, -1, 14, -1, 1, 13, -1,  // Ca40_to_Ar36_He4_approx
        -1, 1, 14, -1, -1, 15, 49,  // Ca40_He4_to_Ti44_approx
        -1, -1, 15, -1, 1, 14, -1,  // Ti44_to_Ca40_He4_approx
        -1, 1, 15, -1, -1, 16, 51,  // Ti44_He4_to_Cr48_approx
        -1, -1, 16, -1, 1, 15, -1,  // Cr48_to_Ti44_He4_approx
        -1, 1, 16, -1, -1, 17, 53,  // Cr48_He4_to_Fe52_approx
        -1, -1, 17, -1, 1, 16, -1,  // Fe52_to_Cr48_He4_approx
        -1, 1, 17, -1, -1, 18, 55,  // Fe52_He4_to_Ni56_approx
        -1, -1, 18, -1, 1, 17, -1,  // Ni56_to_Fe52_He4_approx
        -1, 2, 2, -1, -1, 8, 80,  // C12_C12_to_Mg24_modified
        -1, 4, 4, -1, -1, 12, 81,  // O16_O16_to_S32_modified
        -1, 2, 4, -1, -1, 10, 82,  // C12_O16_to_Si28_modified
        -1, -1, 3, -1, 0, 2, -1,  // N13_to_p_C12_derived
        -1, -1, 4, -1, 1, 2, -1,  // O16_to_He4_C12_derived
        -1, -1, 6, -1, 1, 4, -1,  // Ne20_to_He4_O16_derived
        -1, -1, 8, -1, 0, 7, -1,  // Mg24_to_p_Na23_derived
        -1, -1, 8, -1, 1, 6, -1,  // Mg24_to_He4_Ne20_derived
        -1, -1, 10, -1, 0, 9, -1,  // Si28_to_p_Al27_derived
        -1, -1, 10, -1, 1, 8, -1,  // Si28_to_He4_Mg24_derived
        -1, -1, 12, -1, 0, 11, -1,  // S32_to_p_P31_derived
        -1, -1, 12, -1, 1, 10, -1,  // S32_to_He4_Si28_derived
        -1, -1, 2, 1, 1, 1, -1,  // C12_to_He4_He4_He4_derived
        -1, 0, 4, -1, 1, 3, -1,  // p_O16_to_He4_N13_derived
        -1, 1, 5, -1, 0, 6, 17,  // He4_F17_to_p_Ne20_derived
        -1, 1, 6, -1, 0, 7, -1,  // He4_Ne20_to_p_Na23_derived
        -1, 1, 6, -1, 2, 2, -1,  // He4_Ne20_to_C12_C12_derived
        -1, 0, 7, -1, 2, 2, -1,  // p_Na23_to_C12_C12_derived
        -1, 1, 8, -1, 0, 9, -1,  // He4_Mg24_to_p_Al27_derived
        -1, 1, 8, -1, 2, 4, -1,  // He4_Mg24_to_C12_O16_derived
        -1, 0, 9, -1, 2, 4, -1,  // p_Al27_to_C12_O16_derived
        -1, 1, 10, -1, 0, 11, -1,  // He4_Si28_to_p_P31_derived
        -1, 1, 10, -1, 4, 4, -1,  // He4_Si28_to_O16_O16_derived
        -1, 0, 11, -1, 4, 4, -1,  // p_P31_to_O16_O16_derived
        -1, -1, 8, -1, 2, 2, -1,  // Mg24_to_C12_C12_derived
        -1, -1, 12, -1, 4, 4, -1,  // S32_to_O16_O16_derived
        -1, -1, 10, -1, 2, 4, -1,  // Si28_to_C12_O16_derived
        -1, -1, 5, -1, 0, 4, -1,  // F17_to_O16_p_derived
        -1, -1, -1, -1, -1, -1, -1,  // He4_S32_to_p_Cl35_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ar36_to_He4_S32_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ar36_to_p_Cl35_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Ar36_to_p_K39_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ca40_to_He4_Ar36_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ca40_to_p_K39_removed
        -1, -1, -1, -1, -1, -1, -1,  // He4_Ca40_to_p_Sc43_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ti44_to_He4_Ca40_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ti44_to_p_Sc43_removed
        -1, -1, -1, -1, -1, -1, -1,  // Cr48_to_He4_Ti44_removed
        -1, -1, -1, -1, -1, -1, -1,  // Cr48_to_p_V47_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_V47_to_He4_Ti44_removed
        -1, -1, -1, -1, -1, -1, -1,  // Fe52_to_He4_Cr48_removed
        -1, -1, -1, -1, -1, -1, -1,  // Fe52_to_p_Mn51_removed
        -1, -1, -1, -1, -1, -1, -1,  // p_Mn51_to_He4_Cr48_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ni56_to_He4_Fe52_removed
        -1, -1, -1, -1, -1, -1, -1,  // Ni56_to_p_Co55_removed
        -1, -1, -1, -1, -1, -1, -1  // p_Co55_to_He4_Fe52_removed
    };
}
#endif

void actual_network_init()
{

}
