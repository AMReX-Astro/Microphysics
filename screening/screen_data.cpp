#include <screen_data.H>

#if NUMSCREEN > 0
AMREX_GPU_MANAGED amrex::GpuArray<scrn::screen_factors_t, NSCREEN> scrn::scn_facs;
#endif
