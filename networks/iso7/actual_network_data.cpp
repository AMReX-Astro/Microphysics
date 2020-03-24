#include <actual_network_data.H>

namespace iso7
{
    amrex::GpuArray<amrex::Real,NSpec> bion;
    amrex::GpuArray<amrex::Real,NSpec> mion;
    amrex::GpuArray<amrex::Real,NSpec> wion;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}
