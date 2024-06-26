#include <actual_network.H>


#ifdef NSE_NET
namespace NSE_INDEX
{
    AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 7, Order::C> rate_indices {
        -1, -1, 4, -1, -1, 3, -1,
        -1, -1, 7, -1, -1, 5, -1,
        -1, -1, 8, -1, -1, 6, -1,
        -1, -1, 12, -1, -1, 10, -1,
        -1, -1, 13, -1, -1, 11, -1,
        -1, -1, 15, -1, -1, 13, -1,
        -1, -1, 16, -1, -1, 14, -1,
        -1, -1, 4, -1, 0, 2, -1,
        -1, -1, 5, -1, 0, 3, -1,
        -1, -1, 7, -1, 0, 4, -1,
        -1, -1, 8, -1, 0, 5, -1,
        -1, -1, 9, -1, 0, 6, -1,
        -1, -1, 9, -1, 1, 2, -1,
        -1, -1, 12, -1, 0, 9, -1,
        -1, -1, 13, -1, 0, 10, -1,
        -1, -1, 13, -1, 1, 5, -1,
        -1, -1, 14, -1, 0, 11, -1,
        -1, -1, 14, -1, 1, 6, -1,
        -1, -1, 15, -1, 0, 12, -1,
        -1, -1, 15, -1, 1, 7, -1,
        -1, -1, 16, -1, 0, 13, -1,
        -1, -1, 16, -1, 1, 8, -1,
        -1, -1, 17, -1, 0, 14, -1,
        -1, -1, 17, -1, 1, 9, -1,
        -1, -1, 18, -1, 1, 15, -1,
        -1, -1, 19, -1, 1, 17, -1,
        -1, -1, 2, 1, 1, 1, -1,
        -1, 0, 2, -1, -1, 4, 8,
        -1, 1, 2, -1, -1, 9, 13,
        -1, 0, 3, -1, -1, 5, 9,
        -1, 0, 4, -1, -1, 7, 10,
        -1, 0, 5, -1, -1, 8, 11,
        -1, 1, 5, -1, -1, 13, 16,
        -1, 0, 6, -1, -1, 9, 12,
        -1, 1, 6, -1, -1, 14, 18,
        -1, 1, 7, -1, -1, 15, 20,
        -1, 1, 8, -1, -1, 16, 22,
        -1, 0, 9, -1, -1, 12, 14,
        -1, 1, 9, -1, -1, 17, 24,
        -1, 0, 10, -1, -1, 13, 15,
        -1, 0, 11, -1, -1, 14, 17,
        -1, 0, 12, -1, -1, 15, 19,
        -1, 0, 13, -1, -1, 16, 21,
        -1, 0, 14, -1, -1, 17, 23,
        -1, 1, 15, -1, -1, 18, 25,
        -1, 1, 17, -1, -1, 19, 26,
        -1, 1, 2, -1, 0, 6, -1,
        -1, 2, 2, -1, 1, 17, 65,
        -1, 1, 4, -1, 0, 9, 55,
        -1, 1, 5, -1, 0, 10, -1,
        -1, 0, 6, -1, 1, 2, 47,
        -1, 1, 6, -1, 0, 11, -1,
        -1, 1, 7, -1, 0, 12, 60,
        -1, 1, 8, -1, 0, 13, -1,
        -1, 0, 9, -1, 1, 4, -1,
        -1, 1, 9, -1, 0, 14, -1,
        -1, 2, 9, -1, 1, 19, 66,
        -1, 0, 10, -1, 1, 5, 50,
        -1, 0, 11, -1, 1, 6, 52,
        -1, 0, 12, -1, 1, 7, -1,
        -1, 1, 12, -1, 0, 17, 64,
        -1, 0, 13, -1, 1, 8, 54,
        -1, 0, 14, -1, 1, 9, 56,
        -1, 0, 17, -1, 1, 12, -1,
        -1, 1, 17, -1, 2, 2, -1,
        -1, 1, 19, -1, 2, 9, -1,
        1, 1, 1, -1, -1, 2, 27
    };
}
#endif

void actual_network_init()
{

}
