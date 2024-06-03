#ifndef NSE_TABLE_CHECK_H
#define NSE_TABLE_CHECK_H

#include <iostream>
#include <fstream>
#include <actual_network.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <extern_parameters.H>

#include <burn_type.H>
#include <eos_type.H>


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool in_nse(T& state, const bool relax = false) {
    using namespace Species;

    Real f = relax ? nse_relax_factor : 1.0_rt;
    bool nse_check = false;

    if (state.rho > f * rho_nse && state.T > f * T_nse) {
        Real Fe_group{};
        Real C_group{};
        Real He_group{};
        Real O_group{};
        Real Si_group{};

        auto add_to_group = [&state](Real group, int zion_start, int zion_end) {
            for (int spec = 1; spec <= NumSpec; ++spec) {
                if (NetworkProperties::zion(spec) >= zion_start && NetworkProperties::zion(spec) <= zion_end) {
                    group += state.xn[spec - 1];
                }
            }
        };

        if constexpr (std::is_same_v<T, burn_t>) {
            #ifdef STRANG
            add_to_group(Fe_group, 24, 28); //range for Fe-group?
            add_to_group(C_group, 6, 7);    //range for C-group?
            add_to_group(He_group, 1, 2);   //range for He-group?
            add_to_group(O_group, 8, 8);    //range for O-group?
            add_to_group(Si_group, 14, 14); //range for Si-group?
            #else
            // we need to get the mass fractions from the conserved state passed in
            auto add_to_group_y = [&state](Real group, int zion_start, int zion_end) {
                for (int spec = 1; spec <= NumSpec; ++spec) {
                    if (NetworkProperties::zion(spec) >= zion_start && NetworkProperties::zion(spec) <= zion_end) {
                        group += state.y[SFS + spec - 1];
                    }
                }
                group /= state.rho;
            };
            add_to_group_y(Fe_group, 24, 28); 
            add_to_group_y(C_group, 6, 7);    
            add_to_group_y(He_group, 1, 2);   
            add_to_group_y(O_group, 8, 8);    
            add_to_group_y(Si_group, 14, 14); 
            #endif
        } else {
            // this covers the various eos_t's
            add_to_group(Fe_group, 24, 28); 
            add_to_group(C_group, 6, 7);   
            add_to_group(He_group, 1, 2);   
            add_to_group(O_group, 8, 8);    
            add_to_group(Si_group, 14, 14); 
        }

        if (Fe_group + He_group > f * He_Fe_nse &&
            C_group < C_nse / f &&
            O_group < O_nse / f &&
            Si_group < Si_nse / f) {
            nse_check = true;
        }
    }

    if constexpr (std::is_same<T, burn_t>::value) {
        state.nse = nse_check;
    }

    return nse_check;
}

#endif
