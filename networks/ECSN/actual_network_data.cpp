#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<std::int8_t, 1, Rates::NumRates, 1, 7, amrex::Order::C> rate_indices {
        -1, -1, 5, -1, 1, 2, -1,  // Ne20_to_He4_O16
        -1, 1, 2, -1, -1, 5, 1,  // He4_O16_to_Ne20
        -1, 1, 5, -1, -1, 6, -1,  // He4_Ne20_to_Mg24
        -1, 1, 6, -1, -1, 8, -1,  // He4_Mg24_to_Si28
        -1, 0, 7, -1, -1, 8, -1,  // p_Al27_to_Si28
        -1, 1, 7, -1, -1, 9, -1,  // He4_Al27_to_P31
        -1, 1, 8, -1, -1, 10, -1,  // He4_Si28_to_S32
        -1, 0, 9, -1, -1, 10, -1,  // p_P31_to_S32
        -1, 2, 2, -1, 0, 9, -1,  // O16_O16_to_p_P31
        -1, 2, 2, -1, 1, 8, -1,  // O16_O16_to_He4_Si28
        -1, 1, 6, -1, 0, 7, -1,  // He4_Mg24_to_p_Al27
        -1, 0, 7, -1, 1, 6, 11,  // p_Al27_to_He4_Mg24
        -1, 1, 8, -1, 0, 9, -1,  // He4_Si28_to_p_P31
        -1, 0, 9, -1, 1, 8, 13,  // p_P31_to_He4_Si28
        -1, -1, 4, -1, -1, 5, 17,  // F20_to_Ne20
        -1, -1, 4, -1, -1, 3, -1,  // F20_to_O20
        -1, -1, 5, -1, -1, 4, -1,  // Ne20_to_F20
        -1, -1, 3, -1, -1, 4, 16  // O20_to_F20
    };
}
#endif

void actual_network_init()
{

}
