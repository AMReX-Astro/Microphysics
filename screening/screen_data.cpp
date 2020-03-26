#include <screen_data.H>

using namespace amrex;

namespace scrn {
#if NSCREEN > 0
    AMREX_GPU_MANAGED amrex::GpuArray<screen_factors_t, NSCREEN> scn_facs;
#endif
};
