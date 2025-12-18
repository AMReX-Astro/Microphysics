#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, -1, 8, -1, -1, 7, -1,  // N13_to_C13_reaclib
        -1, -1, 11, -1, -1, 9, -1,  // O14_to_N14_reaclib
        -1, -1, 12, -1, -1, 10, -1,  // O15_to_N15_reaclib
        -1, -1, 15, -1, -1, 14, -1,  // F17_to_O17_reaclib
        -1, -1, 5, -1, 3, 3, -1,  // B8_to_He4_He4_reaclib
        -1, 0, 0, -1, -1, 1, -1,  // p_p_to_d_reaclib_bet_pos
        -1, 0, 0, -1, -1, 1, -1,  // p_p_to_d_reaclib_electron_capture
        -1, 0, 1, -1, -1, 2, -1,  // p_d_to_He3_reaclib
        -1, 1, 1, -1, -1, 3, -1,  // d_d_to_He4_reaclib
        -1, 0, 2, -1, -1, 3, -1,  // p_He3_to_He4_reaclib
        -1, 2, 3, -1, -1, 4, -1,  // He4_He3_to_Be7_reaclib
        -1, 0, 4, -1, -1, 5, -1,  // p_Be7_to_B8_reaclib
        -1, 0, 6, -1, -1, 8, -1,  // p_C12_to_N13_reaclib
        -1, 3, 6, -1, -1, 13, -1,  // He4_C12_to_O16_reaclib
        -1, 0, 7, -1, -1, 9, -1,  // p_C13_to_N14_reaclib
        -1, 0, 8, -1, -1, 11, -1,  // p_N13_to_O14_reaclib
        -1, 0, 9, -1, -1, 12, -1,  // p_N14_to_O15_reaclib
        -1, 3, 9, -1, -1, 16, -1,  // He4_N14_to_F18_reaclib
        -1, 0, 10, -1, -1, 13, -1,  // p_N15_to_O16_reaclib
        -1, 0, 13, -1, -1, 15, -1,  // p_O16_to_F17_reaclib
        -1, 0, 14, -1, -1, 16, -1,  // p_O17_to_F18_reaclib
        -1, 1, 2, -1, 0, 3, -1,  // d_He3_to_p_He4_reaclib
        -1, 3, 8, -1, 0, 13, -1,  // He4_N13_to_p_O16_reaclib
        -1, 0, 10, -1, 3, 6, -1,  // p_N15_to_He4_C12_reaclib
        -1, 3, 11, -1, 0, 15, -1,  // He4_O14_to_p_F17_reaclib
        -1, 0, 14, -1, 3, 9, -1,  // p_O17_to_He4_N14_reaclib
        -1, 0, 16, -1, 3, 12, -1,  // p_F18_to_He4_O15_reaclib
        -1, 2, 2, 0, 0, 3, -1,  // He3_He3_to_p_p_He4_reaclib
        -1, 1, 4, 0, 3, 3, -1,  // d_Be7_to_p_He4_He4_reaclib
        -1, 2, 4, 0, 0, 3, -1,  // He3_Be7_to_p_p_He4_He4_reaclib
        3, 3, 3, -1, -1, 6, -1  // He4_He4_He4_to_C12_reaclib
    };
}
#endif

void actual_network_init()
{

}
