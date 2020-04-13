#include <actual_network_data.H>

namespace iso7
{
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> bion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> mion;
    AMREX_GPU_MANAGED amrex::Array1D<amrex::Real, 1, NumSpec> wion;
}

namespace Rates
{
    amrex::Vector<std::string> names;
}
