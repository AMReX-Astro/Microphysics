#ifndef NEUTRINO_H
#define NEUTRINO_H

// Define available neutrino cooling methods
#define NEUTRINO_METHOD_kipp 1
#define NEUTRINO_METHOD_sneut5 2

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <microphysics_autodiff.H>
#include <fundamental_constants.H>

#include <kipp.H>
#include <sneut5.H>

using namespace amrex::literals;

#if NEUTRINO_METHOD == NEUTRINO_METHOD_kipp
const std::string neutrino_name = "kipp";
#elif NEUTRINO_METHOD == NEUTRINO_METHOD_sneut5
const std::string neutrino_name = "sneut5";
#endif

template <typename number_t>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
number_t neutrino_cooling(const number_t& temp, const amrex::Real dens,
                          const number_t& abar, const number_t& zbar) {

#if NEUTRINO_METHOD == NEUTRINO_METHOD_kipp
    amrex::Real pair, phot, plas, brem;
    return kipp(temp, dens, abar, zbar, pair, phot, plas, brem);

#elif NEUTRINO_METHOD == NEUTRINO_METHOD_sneut5
    amrex::Real pair, phot, plas, brem;
    return sneut5(temp, dens, abar, zbar, pair, phot, plas, brem);
#endif
}

#endif // NEUTRINO_H
