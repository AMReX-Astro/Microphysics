#ifndef NSE_TABLE_SIZE_H
#define NSE_TABLE_SIZE_H

#include <string>

#include <AMReX_REAL.H>

using namespace amrex::literals;

namespace nse_table_size {

    const std::string table_name{"nse_aprox19.tbl"};

    constexpr int ntemp = 101;
    constexpr int nden = 61;
    constexpr int nye = 29;

    constexpr amrex::Real logT_min = 9.4;
    constexpr amrex::Real logT_max = 10.4;
    constexpr amrex::Real dlogT = 0.01;

    constexpr amrex::Real logrho_min = 7.0;
    constexpr amrex::Real logrho_max = 10.0;
    constexpr amrex::Real dlogrho = 0.05;

    constexpr amrex::Real ye_min = 0.43;
    constexpr amrex::Real ye_max = 0.5;
    constexpr amrex::Real dye = 0.0025;

}
#endif
