#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_f20_o20_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 39, 1, 152, 1, 6> j_f20_o20_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 152> j_f20_o20_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 39> j_f20_o20_temp;

    AMREX_GPU_MANAGED table_t j_ne20_f20_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 39, 1, 152, 1, 6> j_ne20_f20_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 152> j_ne20_f20_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 39> j_ne20_f20_temp;

    AMREX_GPU_MANAGED table_t j_o20_f20_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 39, 1, 152, 1, 6> j_o20_f20_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 152> j_o20_f20_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 39> j_o20_f20_temp;

    AMREX_GPU_MANAGED table_t j_f20_ne20_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 39, 1, 152, 1, 6> j_f20_ne20_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 152> j_f20_ne20_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 39> j_f20_ne20_temp;


}


void init_tabular()
{

    std::cout << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_f20_o20_meta.ntemp = 39;
    j_f20_o20_meta.nrhoy = 152;
    j_f20_o20_meta.nvars = 6;
    j_f20_o20_meta.nheader = 5;
    j_f20_o20_meta.file = "20f-20o_electroncapture.dat";

    init_tab_info(j_f20_o20_meta, j_f20_o20_rhoy, j_f20_o20_temp, j_f20_o20_data);


    j_ne20_f20_meta.ntemp = 39;
    j_ne20_f20_meta.nrhoy = 152;
    j_ne20_f20_meta.nvars = 6;
    j_ne20_f20_meta.nheader = 7;
    j_ne20_f20_meta.file = "20ne-20f_electroncapture.dat";

    init_tab_info(j_ne20_f20_meta, j_ne20_f20_rhoy, j_ne20_f20_temp, j_ne20_f20_data);


    j_o20_f20_meta.ntemp = 39;
    j_o20_f20_meta.nrhoy = 152;
    j_o20_f20_meta.nvars = 6;
    j_o20_f20_meta.nheader = 6;
    j_o20_f20_meta.file = "20o-20f_betadecay.dat";

    init_tab_info(j_o20_f20_meta, j_o20_f20_rhoy, j_o20_f20_temp, j_o20_f20_data);


    j_f20_ne20_meta.ntemp = 39;
    j_f20_ne20_meta.nrhoy = 152;
    j_f20_ne20_meta.nvars = 6;
    j_f20_ne20_meta.nheader = 7;
    j_f20_ne20_meta.file = "20f-20ne_betadecay.dat";

    init_tab_info(j_f20_ne20_meta, j_f20_ne20_rhoy, j_f20_ne20_temp, j_f20_ne20_data);



}
