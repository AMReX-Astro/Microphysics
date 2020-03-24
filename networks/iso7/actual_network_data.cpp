#include <actual_network_data.H>

namespace iso7
{
    AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,NSpec> bion;
    AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,NSpec> mion;
    AMREX_GPU_MANAGED amrex::GpuArray<amrex::Real,NSpec> wion;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}
