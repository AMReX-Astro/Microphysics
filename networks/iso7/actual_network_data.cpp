#include <actual_network_data.H>

namespace iso7
{
    AMREX_GPU_MANAGED amrex::Real bion[NSpec];
    AMREX_GPU_MANAGED amrex::Real mion[NSpec];
    AMREX_GPU_MANAGED amrex::Real wion[NSpec];
}

namespace Rates
{
    amrex::Vector<std::string> names;
}
