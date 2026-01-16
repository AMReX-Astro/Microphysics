#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, 1, 3, -1, -1, 4, -1,  // p_C12_to_N13
        -1, 2, 3, -1, -1, 5, -1,  // He4_C12_to_O16
        -1, 2, 5, -1, -1, 6, -1,  // He4_O16_to_Ne20
        -1, 2, 6, -1, -1, 10, -1,  // He4_Ne20_to_Mg24
        -1, 1, 8, -1, -1, 10, -1,  // p_Na23_to_Mg24
        -1, 0, 9, -1, -1, 10, -1,  // n_Mg23_to_Mg24
        -1, 3, 3, -1, 0, 9, -1,  // C12_C12_to_n_Mg23
        -1, 3, 3, -1, 1, 8, -1,  // C12_C12_to_p_Na23
        -1, 3, 3, -1, 2, 6, -1,  // C12_C12_to_He4_Ne20
        -1, 2, 4, -1, 1, 5, -1,  // He4_N13_to_p_O16
        -1, 3, 5, -1, 2, 10, -1,  // C12_O16_to_He4_Mg24
        -1, 2, 6, -1, 1, 8, -1,  // He4_Ne20_to_p_Na23
        -1, 1, 8, -1, 2, 6, 12,  // p_Na23_to_He4_Ne20
        -1, 0, 9, -1, 1, 8, -1,  // n_Mg23_to_p_Na23
        -1, 0, 9, -1, 2, 6, -1,  // n_Mg23_to_He4_Ne20
        -1, 0, 9, -1, 3, 3, 7,  // n_Mg23_to_C12_C12
        2, 2, 2, -1, -1, 3, -1,  // He4_He4_He4_to_C12
        -1, -1, 8, -1, -1, 7, -1,  // Na23_to_Ne23
        -1, -1, 7, -1, -1, 8, 18,  // Ne23_to_Na23
        -1, -1, 9, -1, -1, 8, -1,  // Mg23_to_Na23
        -1, -1, 0, -1, -1, 1, 22,  // n_to_p
        -1, -1, 1, -1, -1, 0, -1  // p_to_n
    };
}
#endif

void actual_network_init()
{

}
