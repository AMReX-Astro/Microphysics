#include <actual_network_data.H>

namespace iso7
{
    AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,NumSpec> bion;
    AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,NumSpec> mion;
    AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,NumSpec> wion;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}
