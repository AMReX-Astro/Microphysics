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
        -1, -1, 21, -1, -1, 19, -1,
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
        -1, -1, 18, -1, 1, 10, -1,
        -1, -1, 19, -1, 0, 18, -1,
        -1, -1, 19, -1, 1, 13, -1,
        -1, -1, 20, -1, 1, 14, -1,
        -1, -1, 21, -1, 1, 15, -1,
        -1, -1, 22, -1, 0, 20, -1,
        -1, -1, 22, -1, 1, 17, -1,
        -1, -1, 24, -1, 0, 23, -1,
        -1, -1, 24, -1, 1, 22, -1,
        -1, -1, 26, -1, 0, 25, -1,
        -1, -1, 26, -1, 1, 24, -1,
        -1, -1, 2, 1, 1, 1, -1,
        -1, 0, 2, -1, -1, 4, 9,
        -1, 1, 2, -1, -1, 9, 14,
        -1, 0, 3, -1, -1, 5, 10,
        -1, 0, 4, -1, -1, 7, 11,
        -1, 0, 5, -1, -1, 8, 12,
        -1, 1, 5, -1, -1, 13, 17,
        -1, 0, 6, -1, -1, 9, 13,
        -1, 1, 6, -1, -1, 14, 19,
        -1, 1, 7, -1, -1, 15, 21,
        -1, 1, 8, -1, -1, 16, 23,
        -1, 0, 9, -1, -1, 12, 15,
        -1, 1, 9, -1, -1, 17, 25,
        -1, 0, 10, -1, -1, 13, 16,
        -1, 1, 10, -1, -1, 18, 26,
        -1, 0, 11, -1, -1, 14, 18,
        -1, 0, 12, -1, -1, 15, 20,
        -1, 0, 13, -1, -1, 16, 22,
        -1, 1, 13, -1, -1, 19, 28,
        -1, 0, 14, -1, -1, 17, 24,
        -1, 1, 14, -1, -1, 20, 29,
        -1, 1, 15, -1, -1, 21, 30,
        -1, 1, 17, -1, -1, 22, 32,
        -1, 0, 18, -1, -1, 19, 27,
        -1, 0, 20, -1, -1, 22, 31,
        -1, 1, 22, -1, -1, 24, 34,
        -1, 0, 23, -1, -1, 24, 33,
        -1, 1, 24, -1, -1, 26, 36,
        -1, 0, 25, -1, -1, 26, 35,
        -1, 1, 2, -1, 0, 6, -1,
        -1, 2, 2, -1, 0, 20, -1,
        -1, 2, 2, -1, 1, 17, -1,
        -1, 1, 4, -1, 0, 9, 75,
        -1, 1, 5, -1, 0, 10, -1,
        -1, 0, 6, -1, 1, 2, 66,
        -1, 1, 6, -1, 0, 11, -1,
        -1, 1, 7, -1, 0, 12, 83,
        -1, 1, 8, -1, 0, 13, -1,
        -1, 0, 9, -1, 1, 4, -1,
        -1, 1, 9, -1, 0, 14, -1,
        -1, 2, 9, -1, 0, 23, -1,
        -1, 2, 9, -1, 1, 22, -1,
        -1, 9, 9, -1, 0, 25, -1,
        -1, 9, 9, -1, 1, 24, -1,
        -1, 0, 10, -1, 1, 5, 70,
        -1, 0, 11, -1, 1, 6, 72,
        -1, 0, 12, -1, 1, 7, -1,
        -1, 1, 12, -1, 0, 17, 89,
        -1, 0, 13, -1, 1, 8, 74,
        -1, 1, 13, -1, 0, 18, 91,
        -1, 0, 14, -1, 1, 9, 76,
        -1, 1, 16, -1, 0, 19, 92,
        -1, 0, 17, -1, 1, 12, -1,
        -1, 1, 17, -1, 0, 20, -1,
        -1, 0, 18, -1, 1, 13, -1,
        -1, 0, 19, -1, 1, 16, -1,
        -1, 0, 20, -1, 1, 17, 90,
        -1, 1, 22, -1, 0, 23, -1,
        -1, 0, 23, -1, 1, 22, 94,
        -1, 1, 24, -1, 0, 25, -1,
        -1, 0, 25, -1, 1, 24, 96,
        1, 1, 1, -1, -1, 2, 37,
        -1, 2, 2, -1, -1, 22, -1,
        -1, 9, 9, -1, -1, 26, -1,
        -1, 2, 9, -1, -1, 24, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1,
        -1, 1, 26, -1, -1, 27, 139,
        -1, -1, 27, -1, 1, 26, -1,
        -1, 1, 27, -1, -1, 28, 141,
        -1, -1, 28, -1, 1, 27, -1,
        -1, 1, 28, -1, -1, 29, 143,
        -1, -1, 29, -1, 1, 28, -1,
        -1, 1, 29, -1, -1, 30, 145,
        -1, -1, 30, -1, 1, 29, -1,
        -1, 1, 30, -1, -1, 31, 147,
        -1, -1, 31, -1, 1, 30, -1,
        -1, 1, 31, -1, -1, 32, 149,
        -1, -1, 32, -1, 1, 31, -1
    };
}
#endif

void actual_network_init()
{

}
