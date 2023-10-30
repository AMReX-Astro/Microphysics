#ifndef NSE_TABLE_SIZE_H
#define NSE_TABLE_SIZE_H

#include <string>

namespace nse_table_size {

    const std::string table_name{"nse_aprox19.tbl"};

    constexpr int ntemp = 51;
    constexpr int nden = 31;
    constexpr int nye = 15;

    constexpr Real logT_min = 9.4;
    constexpr Real logT_max = 10.4;
    constexpr Real dlogT = 0.02;

    constexpr Real logrho_min = 7.0;
    constexpr Real logrho_max = 10.0;
    constexpr Real dlogrho = 0.1;

    constexpr Real ye_min = 0.43;
    constexpr Real ye_max = 0.5;
    constexpr Real dye = 0.005;

}
#endif
