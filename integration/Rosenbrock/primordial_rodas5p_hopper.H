#ifndef PRIMORDIAL_RODAS5P_HOPPER_H
#define PRIMORDIAL_RODAS5P_HOPPER_H

#include <AMReX_REAL.H>

#include <burn_type.H>

#ifdef AMREX_USE_CUDA

namespace primordial_rodas5p_hopper {

constexpr int zones_per_block = 16;
constexpr int warps_per_block = 8;
constexpr int threads_per_block = warps_per_block * 32;
constexpr int neq = neqs;
constexpr int scratch_count = 534;

constexpr std::size_t align_up (const std::size_t value, const std::size_t alignment)
{
    return ((value + alignment - 1) / alignment) * alignment;
}

constexpr std::size_t bytes_for_reals (const std::size_t n)
{
    return align_up(n * sizeof(amrex::Real), 128);
}

constexpr std::size_t bytes_for_ints (const std::size_t n)
{
    return align_up(n * sizeof(int), 128);
}

constexpr std::size_t shared_bytes ()
{
    std::size_t bytes = 0;
    bytes += bytes_for_reals(neq * zones_per_block);                 // y
    bytes += bytes_for_reals(neq * zones_per_block);                 // ynew
    bytes += bytes_for_reals(neq * zones_per_block);                 // work
    bytes += bytes_for_reals(8 * neq * zones_per_block);             // ak
    bytes += bytes_for_reals(neq * neq * zones_per_block);           // jac
    bytes += bytes_for_reals(neq * zones_per_block);                 // rhs_tmp
    bytes += bytes_for_reals(scratch_count * zones_per_block);       // generated expression scratch
    bytes += bytes_for_ints(neq * zones_per_block);                  // pivots
    bytes += bytes_for_ints(4 * zones_per_block);                    // status/counters
    return bytes;
}

void burn_device (burn_t* states, amrex::Real dt, int nstates);

} // namespace primordial_rodas5p_hopper

#endif

#endif
