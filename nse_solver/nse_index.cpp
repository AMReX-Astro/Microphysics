#include <nse_index.H>

// Declare extern nse index variables

bool NSE_INDEX::initialized;
AMREX_GPU_MANAGED int NSE_INDEX::p_index {-1};
AMREX_GPU_MANAGED int NSE_INDEX::h1_index {-1};
AMREX_GPU_MANAGED int NSE_INDEX::n_index {-1};
AMREX_GPU_MANAGED int NSE_INDEX::he4_index {-1};
AMREX_GPU_MANAGED amrex::Array2D<int, 1, Rates::NumRates, 1, 6> NSE_INDEX::rate_indices;
