#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <actual_network.H>

using namespace reaclib_rates;
using namespace Rates;

struct tf_t {
    Real T9;
    Real T9i;
    Real T943i;
    Real T923i;
    Real T913i;
    Real T913;
    Real T923;
    Real T953;
    Real lnT9;
};

AMREX_INLINE
void init_reaclib() {

    const std::string secret_code_reference = "SZG2JEFM4E0X9HDIBZDL4A9ONAKI24WF";

    amrex::Print() << "reading in reaclib rates..." << std::endl;

    std::ifstream table;
    table.open("reaclib_rate_metadata.dat");

    std::string secret_code_file;
    std::getline(table, secret_code_file);

    if (secret_code_file != secret_code_reference) {
        amrex::Error("invalid reaclib_rate_metadata.dat");
    }

    std::string line;

    // read in all of the reaclib rate sets

    for (int ireaclib = 1; ireaclib <= NumReaclibSets; ++ireaclib) {
        for (int icoeff = 1; icoeff <= 7; ++icoeff) {
            std::getline(table, line);
            std::istringstream data(line);
            data >> ctemp_rate(icoeff, ireaclib);
        }
    }

    for (int ireaclib = 1; ireaclib <= NrateReaclib; ++ireaclib) {
        std::getline(table, line);
        std::istringstream data(line);
        data >> rate_start_idx(ireaclib);
    }

    for (int ireaclib = 1; ireaclib <= NrateReaclib; ++ireaclib) {
        std::getline(table, line);
        std::istringstream data(line);
        data >> rate_extra_mult(ireaclib);
    }

    table.close();

}

AMREX_INLINE
void net_screening_init() {

    using namespace Species;

    // note: we need to set these up in the same order that we evaluate the
    // rates in actual_rhs.H (yes, it's ugly)

    int jscr = 0;

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[H1-1], aion[H1-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[H2-1], aion[H2-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[H2-1], aion[H2-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[He3-1], aion[He3-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He3-1], aion[He3-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Li6-1], aion[Li6-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Li6-1], aion[Li6-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Li7-1], aion[Li7-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Be7-1], aion[Be7-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Be9-1], aion[Be9-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[B11-1], aion[B11-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[C12-1], aion[C12-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[C12-1], aion[C12-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[C13-1], aion[C13-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[C14-1], aion[C14-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[C14-1], aion[C14-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[N13-1], aion[N13-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[N14-1], aion[N14-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[N14-1], aion[N14-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[N15-1], aion[N15-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[N15-1], aion[N15-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[O14-1], aion[O14-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[O15-1], aion[O15-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[O17-1], aion[O17-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[O17-1], aion[O17-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[O18-1], aion[O18-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[O18-1], aion[O18-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[F17-1], aion[F17-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[F17-1], aion[F17-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[F18-1], aion[F18-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[F18-1], aion[F18-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[F19-1], aion[F19-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[F19-1], aion[F19-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ne19-1], aion[Ne19-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ne20-1], aion[Ne20-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ne20-1], aion[Ne20-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ne21-1], aion[Ne21-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ne21-1], aion[Ne21-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ne22-1], aion[Ne22-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ne22-1], aion[Ne22-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Na21-1], aion[Na21-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Na22-1], aion[Na22-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Na22-1], aion[Na22-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Na23-1], aion[Na23-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Na23-1], aion[Na23-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mg24-1], aion[Mg24-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mg24-1], aion[Mg24-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mg25-1], aion[Mg25-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mg25-1], aion[Mg25-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mg26-1], aion[Mg26-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mg26-1], aion[Mg26-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Al25-1], aion[Al25-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Al26-1], aion[Al26-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Al27-1], aion[Al27-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Al27-1], aion[Al27-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Si28-1], aion[Si28-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Si28-1], aion[Si28-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Si29-1], aion[Si29-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Si29-1], aion[Si29-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Si30-1], aion[Si30-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Si30-1], aion[Si30-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Si31-1], aion[Si31-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Si31-1], aion[Si31-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Si32-1], aion[Si32-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Si32-1], aion[Si32-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[P29-1], aion[P29-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[P30-1], aion[P30-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[P31-1], aion[P31-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[P31-1], aion[P31-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[P32-1], aion[P32-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[P32-1], aion[P32-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[P33-1], aion[P33-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[P33-1], aion[P33-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[S32-1], aion[S32-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[S32-1], aion[S32-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[S33-1], aion[S33-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[S33-1], aion[S33-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[S34-1], aion[S34-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[S34-1], aion[S34-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[S35-1], aion[S35-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[S35-1], aion[S35-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[S36-1], aion[S36-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[S36-1], aion[S36-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cl33-1], aion[Cl33-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cl34-1], aion[Cl34-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cl35-1], aion[Cl35-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cl35-1], aion[Cl35-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cl36-1], aion[Cl36-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cl36-1], aion[Cl36-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cl37-1], aion[Cl37-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cl37-1], aion[Cl37-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ar36-1], aion[Ar36-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ar36-1], aion[Ar36-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ar37-1], aion[Ar37-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ar37-1], aion[Ar37-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ar38-1], aion[Ar38-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ar38-1], aion[Ar38-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ar39-1], aion[Ar39-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ar39-1], aion[Ar39-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ar40-1], aion[Ar40-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ar40-1], aion[Ar40-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[K39-1], aion[K39-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[K39-1], aion[K39-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[K40-1], aion[K40-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[K40-1], aion[K40-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[K41-1], aion[K41-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[K41-1], aion[K41-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca40-1], aion[Ca40-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca41-1], aion[Ca41-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca42-1], aion[Ca42-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca42-1], aion[Ca42-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca43-1], aion[Ca43-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca43-1], aion[Ca43-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca44-1], aion[Ca44-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca44-1], aion[Ca44-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca45-1], aion[Ca45-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca45-1], aion[Ca45-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca46-1], aion[Ca46-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca46-1], aion[Ca46-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca47-1], aion[Ca47-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca47-1], aion[Ca47-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca48-1], aion[Ca48-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Sc43-1], aion[Sc43-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Sc43-1], aion[Sc43-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Sc44-1], aion[Sc44-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Sc44-1], aion[Sc44-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Sc45-1], aion[Sc45-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Sc45-1], aion[Sc45-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Sc46-1], aion[Sc46-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Sc46-1], aion[Sc46-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Sc47-1], aion[Sc47-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Sc47-1], aion[Sc47-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Sc48-1], aion[Sc48-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Sc48-1], aion[Sc48-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Sc49-1], aion[Sc49-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti44-1], aion[Ti44-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ti45-1], aion[Ti45-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti45-1], aion[Ti45-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ti46-1], aion[Ti46-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti46-1], aion[Ti46-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ti47-1], aion[Ti47-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti47-1], aion[Ti47-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ti48-1], aion[Ti48-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti48-1], aion[Ti48-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ti49-1], aion[Ti49-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti49-1], aion[Ti49-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ti50-1], aion[Ti50-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti50-1], aion[Ti50-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ti51-1], aion[Ti51-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[V46-1], aion[V46-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[V47-1], aion[V47-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[V47-1], aion[V47-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[V48-1], aion[V48-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[V48-1], aion[V48-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[V49-1], aion[V49-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[V49-1], aion[V49-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[V50-1], aion[V50-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[V50-1], aion[V50-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[V51-1], aion[V51-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[V51-1], aion[V51-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[V52-1], aion[V52-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cr48-1], aion[Cr48-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cr49-1], aion[Cr49-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cr49-1], aion[Cr49-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cr50-1], aion[Cr50-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cr50-1], aion[Cr50-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cr51-1], aion[Cr51-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cr51-1], aion[Cr51-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cr52-1], aion[Cr52-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cr52-1], aion[Cr52-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cr53-1], aion[Cr53-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cr53-1], aion[Cr53-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cr54-1], aion[Cr54-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cr54-1], aion[Cr54-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mn50-1], aion[Mn50-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mn51-1], aion[Mn51-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mn51-1], aion[Mn51-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mn52-1], aion[Mn52-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mn52-1], aion[Mn52-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mn53-1], aion[Mn53-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mn53-1], aion[Mn53-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mn54-1], aion[Mn54-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mn54-1], aion[Mn54-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Mn55-1], aion[Mn55-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mn55-1], aion[Mn55-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Fe52-1], aion[Fe52-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Fe52-1], aion[Fe52-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Fe53-1], aion[Fe53-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Fe53-1], aion[Fe53-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Fe54-1], aion[Fe54-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Fe54-1], aion[Fe54-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Fe55-1], aion[Fe55-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Fe55-1], aion[Fe55-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Fe56-1], aion[Fe56-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Fe56-1], aion[Fe56-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Fe57-1], aion[Fe57-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Fe57-1], aion[Fe57-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Fe58-1], aion[Fe58-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Fe58-1], aion[Fe58-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Co53-1], aion[Co53-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Co54-1], aion[Co54-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Co55-1], aion[Co55-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Co55-1], aion[Co55-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Co56-1], aion[Co56-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Co56-1], aion[Co56-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Co57-1], aion[Co57-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Co57-1], aion[Co57-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Co58-1], aion[Co58-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Co58-1], aion[Co58-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Co59-1], aion[Co59-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Co59-1], aion[Co59-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni56-1], aion[Ni56-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni56-1], aion[Ni56-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni57-1], aion[Ni57-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni57-1], aion[Ni57-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni58-1], aion[Ni58-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni58-1], aion[Ni58-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni59-1], aion[Ni59-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni59-1], aion[Ni59-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni60-1], aion[Ni60-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni60-1], aion[Ni60-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni61-1], aion[Ni61-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni61-1], aion[Ni61-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni62-1], aion[Ni62-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni62-1], aion[Ni62-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni63-1], aion[Ni63-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ni64-1], aion[Ni64-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu58-1], aion[Cu58-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cu58-1], aion[Cu58-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu59-1], aion[Cu59-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cu59-1], aion[Cu59-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu60-1], aion[Cu60-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cu60-1], aion[Cu60-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu61-1], aion[Cu61-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu62-1], aion[Cu62-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu63-1], aion[Cu63-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu64-1], aion[Cu64-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Cu65-1], aion[Cu65-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Zn59-1], aion[Zn59-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Zn60-1], aion[Zn60-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Zn61-1], aion[Zn61-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Zn62-1], aion[Zn62-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Zn63-1], aion[Zn63-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ga62-1], aion[Ga62-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ga63-1], aion[Ga63-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[He3-1], aion[He3-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[Li6-1], aion[Li6-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Li7-1], aion[Li7-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Be7-1], aion[Be7-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Be9-1], aion[Be9-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[B10-1], aion[B10-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[B10-1], aion[B10-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[B11-1], aion[B11-1]);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[C12-1], aion[C12-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[C13-1], aion[C13-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[C13-1], aion[C13-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[C14-1], aion[C14-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[N13-1], aion[N13-1]);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[O16-1], aion[O16-1], zion[O16-1], aion[O16-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ne18-1], aion[Ne18-1]);

    add_screening_factor(jscr++, zion[C12-1], aion[C12-1], zion[Ne20-1], aion[Ne20-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Na21-1], aion[Na21-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Mg23-1], aion[Mg23-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Al26-1], aion[Al26-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[K37-1], aion[K37-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[K38-1], aion[K38-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca40-1], aion[Ca40-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ca41-1], aion[Ca41-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ca48-1], aion[Ca48-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Sc49-1], aion[Sc49-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ti51-1], aion[Ti51-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[V52-1], aion[V52-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Ni63-1], aion[Ni63-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cu57-1], aion[Cu57-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cu61-1], aion[Cu61-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cu62-1], aion[Cu62-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Cu63-1], aion[Cu63-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Zn60-1], aion[Zn60-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[Zn61-1], aion[Zn61-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Zn64-1], aion[Zn64-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Zn65-1], aion[Zn65-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Zn66-1], aion[Zn66-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[Ga64-1], aion[Ga64-1]);

    add_screening_factor(jscr++, zion[He3-1], aion[He3-1], zion[He3-1], aion[He3-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[Li7-1], aion[Li7-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[Be7-1], aion[Be7-1]);

    add_screening_factor(jscr++, zion[He3-1], aion[He3-1], zion[Li7-1], aion[Li7-1]);

    add_screening_factor(jscr++, zion[He3-1], aion[He3-1], zion[Be7-1], aion[Be7-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[He4-1], aion[He4-1], 4.0_rt, 8.0_rt);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[H1-1], aion[H1-1]);

    add_screening_factor(jscr++, zion[H1-1], aion[H1-1], zion[He4-1], aion[He4-1]);

    add_screening_factor(jscr++, zion[H2-1], aion[H2-1], zion[He4-1], aion[He4-1]);


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
tf_t evaluate_tfactors(const Real T)
{

    tf_t tf;
    tf.T9 = T / 1.e9_rt;
    tf.T9i = 1.0_rt / tf.T9;
    tf.T913i = std::pow(tf.T9i, 1.0_rt/3.0_rt);
    tf.T923i = tf.T913i * tf.T913i;
    tf.T943i = tf.T9i * tf.T913i;
    tf.T913 = std::pow(tf.T9, 1.0_rt/3.0_rt);
    tf.T923 = tf.T913 * tf.T913;
    tf.T953 = tf.T9 * tf.T923;
    tf.lnT9 = std::log(tf.T9);

    return tf;
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void reaclib_evaluate(plasma_state_t& pstate, const tf_t& tfactors, const int iwhich,
                      Real& rate, Real& drate_dt) {

    // iwhich is the rate to evaluate
    // we then return the rate and its temperature derivative

    Real irate = 0.0e0_rt;
    Real dirate_dt = 0.0e0_rt;

    rate = 0.0e0_rt;
    drate_dt = 0.0e0_rt;

    // Get the number of additional Reaclib sets for this rate
    // Total number of Reaclib sets for this rate is m + 1
    int m = reaclib_rates::rate_extra_mult(iwhich);

    int istart = reaclib_rates::rate_start_idx(iwhich);

    for (int i = 0; i <= m; ++i) {

        // compute the log of the i component of the rate, in Reaclib form

        Real lnirate =
            reaclib_rates::ctemp_rate(1, istart+i) +
            reaclib_rates::ctemp_rate(2, istart+i) * tfactors.T9i +
            reaclib_rates::ctemp_rate(3, istart+i) * tfactors.T913i +
            reaclib_rates::ctemp_rate(4, istart+i) * tfactors.T913 +
            reaclib_rates::ctemp_rate(5, istart+i) * tfactors.T9 +
            reaclib_rates::ctemp_rate(6, istart+i) * tfactors.T953 +
            reaclib_rates::ctemp_rate(7, istart+i) * tfactors.lnT9;

        // compute the derivative of the log of the rate

        Real dlnirate_dT9{0.0};

        if constexpr (do_T_derivatives == 1) {
            dlnirate_dT9 =
                -reaclib_rates::ctemp_rate(2, istart+i) * tfactors.T9i * tfactors.T9i +
                -(1.0_rt/3.0_rt) * reaclib_rates::ctemp_rate(3, istart+i) * tfactors.T943i +
                (1.0_rt/3.0_rt) * reaclib_rates::ctemp_rate(4, istart+i) * tfactors.T923i +
                reaclib_rates::ctemp_rate(5, istart+i) +
                (5.0_rt/3.0_rt) * reaclib_rates::ctemp_rate(6, istart+i) * tfactors.T923 +
                reaclib_rates::ctemp_rate(7, istart+i) * tfactors.T9i;
        }

        // If the rate will be in the approx. interval [0.0, 1.0E-100], replace by 0.0
        // This avoids issues with passing very large negative values to EXP
        // and getting results between 0.0 and 1.0E-308, the limit for IEEE 754.
        // And avoids SIGFPE in CVODE due to tiny rates.

        lnirate = amrex::max(lnirate, -230.0e0_rt);
        irate = std::exp(lnirate);
        rate += irate;

        if constexpr (do_T_derivatives == 1) {
            dirate_dt = irate * dlnirate_dT9 / 1.0e9_rt;
            drate_dt += dirate_dt;
        }
    }
}

#endif
