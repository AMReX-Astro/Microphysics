#ifndef NSE_TABLE_TYPE_H
#define NSE_TABLE_TYPE_H

#include <AMReX_REAL.H>
#include <network.H>

struct nse_table_t {

    amrex::Real T;
    amrex::Real rho;
    amrex::Real Ye;
    amrex::Real abar;
    amrex::Real bea;
    amrex::Real dyedt;
    amrex::Real dabardt;
    amrex::Real dbeadt;
    amrex::Real e_nu;
    amrex::Real X[NumSpec];

};


#endif
