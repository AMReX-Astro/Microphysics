#include <AMReX_Array.H>
#include <string>
#include <table_rates.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace rate_tables
{

    AMREX_GPU_MANAGED table_t j_co55_fe55_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_co55_fe55_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_co55_fe55_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_co55_fe55_temp;

    AMREX_GPU_MANAGED table_t j_co56_fe56_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_co56_fe56_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_co56_fe56_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_co56_fe56_temp;

    AMREX_GPU_MANAGED table_t j_co56_ni56_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_co56_ni56_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_co56_ni56_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_co56_ni56_temp;

    AMREX_GPU_MANAGED table_t j_co57_ni57_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_co57_ni57_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_co57_ni57_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_co57_ni57_temp;

    AMREX_GPU_MANAGED table_t j_fe55_co55_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_fe55_co55_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_fe55_co55_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_fe55_co55_temp;

    AMREX_GPU_MANAGED table_t j_fe55_mn55_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_fe55_mn55_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_fe55_mn55_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_fe55_mn55_temp;

    AMREX_GPU_MANAGED table_t j_fe56_co56_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_fe56_co56_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_fe56_co56_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_fe56_co56_temp;

    AMREX_GPU_MANAGED table_t j_mn55_fe55_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_mn55_fe55_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_mn55_fe55_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_mn55_fe55_temp;

    AMREX_GPU_MANAGED table_t j_ni56_co56_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_ni56_co56_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_ni56_co56_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_ni56_co56_temp;

    AMREX_GPU_MANAGED table_t j_ni57_co57_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_ni57_co57_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_ni57_co57_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_ni57_co57_temp;

    AMREX_GPU_MANAGED table_t j_n_p_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_n_p_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_n_p_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_n_p_temp;

    AMREX_GPU_MANAGED table_t j_p_n_meta;
    AMREX_GPU_MANAGED Array3D<Real, 1, 13, 1, 11, 1, 6> j_p_n_data;
    AMREX_GPU_MANAGED Array1D<Real, 1, 11> j_p_n_rhoy;
    AMREX_GPU_MANAGED Array1D<Real, 1, 13> j_p_n_temp;


}


void init_tabular()
{

    amrex::Print() << "reading in network electron-capture / beta-decay tables..." << std::endl;

    using namespace rate_tables;

    j_co55_fe55_meta.ntemp = 13;
    j_co55_fe55_meta.nrhoy = 11;
    j_co55_fe55_meta.nvars = 6;
    j_co55_fe55_meta.nheader = 5;

    init_tab_info(j_co55_fe55_meta, "55co-55fe_electroncapture.dat", j_co55_fe55_rhoy, j_co55_fe55_temp, j_co55_fe55_data);


    j_co56_fe56_meta.ntemp = 13;
    j_co56_fe56_meta.nrhoy = 11;
    j_co56_fe56_meta.nvars = 6;
    j_co56_fe56_meta.nheader = 5;

    init_tab_info(j_co56_fe56_meta, "56co-56fe_electroncapture.dat", j_co56_fe56_rhoy, j_co56_fe56_temp, j_co56_fe56_data);


    j_co56_ni56_meta.ntemp = 13;
    j_co56_ni56_meta.nrhoy = 11;
    j_co56_ni56_meta.nvars = 6;
    j_co56_ni56_meta.nheader = 5;

    init_tab_info(j_co56_ni56_meta, "56co-56ni_betadecay.dat", j_co56_ni56_rhoy, j_co56_ni56_temp, j_co56_ni56_data);


    j_co57_ni57_meta.ntemp = 13;
    j_co57_ni57_meta.nrhoy = 11;
    j_co57_ni57_meta.nvars = 6;
    j_co57_ni57_meta.nheader = 5;

    init_tab_info(j_co57_ni57_meta, "57co-57ni_betadecay.dat", j_co57_ni57_rhoy, j_co57_ni57_temp, j_co57_ni57_data);


    j_fe55_co55_meta.ntemp = 13;
    j_fe55_co55_meta.nrhoy = 11;
    j_fe55_co55_meta.nvars = 6;
    j_fe55_co55_meta.nheader = 5;

    init_tab_info(j_fe55_co55_meta, "55fe-55co_betadecay.dat", j_fe55_co55_rhoy, j_fe55_co55_temp, j_fe55_co55_data);


    j_fe55_mn55_meta.ntemp = 13;
    j_fe55_mn55_meta.nrhoy = 11;
    j_fe55_mn55_meta.nvars = 6;
    j_fe55_mn55_meta.nheader = 5;

    init_tab_info(j_fe55_mn55_meta, "55fe-55mn_electroncapture.dat", j_fe55_mn55_rhoy, j_fe55_mn55_temp, j_fe55_mn55_data);


    j_fe56_co56_meta.ntemp = 13;
    j_fe56_co56_meta.nrhoy = 11;
    j_fe56_co56_meta.nvars = 6;
    j_fe56_co56_meta.nheader = 5;

    init_tab_info(j_fe56_co56_meta, "56fe-56co_betadecay.dat", j_fe56_co56_rhoy, j_fe56_co56_temp, j_fe56_co56_data);


    j_mn55_fe55_meta.ntemp = 13;
    j_mn55_fe55_meta.nrhoy = 11;
    j_mn55_fe55_meta.nvars = 6;
    j_mn55_fe55_meta.nheader = 5;

    init_tab_info(j_mn55_fe55_meta, "55mn-55fe_betadecay.dat", j_mn55_fe55_rhoy, j_mn55_fe55_temp, j_mn55_fe55_data);


    j_ni56_co56_meta.ntemp = 13;
    j_ni56_co56_meta.nrhoy = 11;
    j_ni56_co56_meta.nvars = 6;
    j_ni56_co56_meta.nheader = 5;

    init_tab_info(j_ni56_co56_meta, "56ni-56co_electroncapture.dat", j_ni56_co56_rhoy, j_ni56_co56_temp, j_ni56_co56_data);


    j_ni57_co57_meta.ntemp = 13;
    j_ni57_co57_meta.nrhoy = 11;
    j_ni57_co57_meta.nvars = 6;
    j_ni57_co57_meta.nheader = 5;

    init_tab_info(j_ni57_co57_meta, "57ni-57co_electroncapture.dat", j_ni57_co57_rhoy, j_ni57_co57_temp, j_ni57_co57_data);


    j_n_p_meta.ntemp = 13;
    j_n_p_meta.nrhoy = 11;
    j_n_p_meta.nvars = 6;
    j_n_p_meta.nheader = 5;

    init_tab_info(j_n_p_meta, "n-p_betadecay.dat", j_n_p_rhoy, j_n_p_temp, j_n_p_data);


    j_p_n_meta.ntemp = 13;
    j_p_n_meta.nrhoy = 11;
    j_p_n_meta.nvars = 6;
    j_p_n_meta.nheader = 5;

    init_tab_info(j_p_n_meta, "p-n_electroncapture.dat", j_p_n_rhoy, j_p_n_temp, j_p_n_data);



}
