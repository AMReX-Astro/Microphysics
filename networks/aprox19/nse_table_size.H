#ifndef NSE_TABLE_SIZE_H
#define NSE_TABLE_SIZE_H

#include <string>

#include <AMReX_REAL.H>

using namespace amrex;

namespace nse_table_size {

    const std::string table_name{"nse_aprox19.tbl"};

    constexpr int ntemp = 101;
    constexpr int nden = 61;
    constexpr int nye = 29;

    constexpr Real logT_min = 9.4;
    constexpr Real logT_max = 10.4;
    constexpr Real dlogT = 0.01;

    constexpr Real logrho_min = 7.0;
    constexpr Real logrho_max = 10.0;
    constexpr Real dlogrho = 0.05;

    constexpr Real ye_min = 0.43;
    constexpr Real ye_max = 0.5;
    constexpr Real dye = 0.0025;

}
#endif
