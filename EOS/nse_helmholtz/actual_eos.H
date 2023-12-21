#ifndef ACTUAL_EOS_H
#define ACTUAL_EOS_H

#include <iostream>
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_ParallelDescriptor.H>
#include <extern_parameters.H>
#include <eos_type.H>
#include <core_eos_interface.H>
#include <nse_table_check.H>
#include <nse_eos.H>

template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_eos (I input, T& state)
{

    if (in_nse(state) && (input == eos_input_re || input == eos_input_rp)) {


        if (input == eos_input_rp) {

            if constexpr (has_pressure<T>::value) {
                nse_T_abar_from_p(state.rho, state.p, state.aux[iye],
                                  state.T, state.aux[iabar]);
            }

        } else {

            nse_T_abar_from_e(state.rho, state.e, state.aux[iye],
                              state.T, state.aux[iabar]);

        }

        // this just solved for T / Abar, now get the rest of the state
        core_eos_interface(eos_input_rt, state);

    } else {
        core_eos_interface(input, state);
    }

}


AMREX_INLINE
void actual_eos_init ()
{
    using namespace helmholtz;

    Real dth, dt2, dti, dt2i;
    Real dd, dd2, ddi, dd2i;

    // Read in the runtime parameters

    input_is_constant = eos_input_is_constant;
    do_coulomb = use_eos_coulomb;
    ttol = eos_ttol;
    dtol = eos_dtol;

    //    read the helmholtz free energy table
    tlo   = 3.0e0_rt;
    thi   = 13.0e0_rt;
    tstp  = (thi - tlo) / ((Real) (jmax-1));
    tstpi = 1.0e0_rt / tstp;
    dlo   = -12.0e0_rt;
    dhi   = 15.0e0_rt;
    dstp  = (dhi - dlo) / ((Real) (imax-1));
    dstpi = 1.0e0_rt / dstp;

    for (int j = 0; j < jmax; ++j) {
        Real tsav = tlo + j * tstp;
        t[j] = std::pow(10.0e0_rt, tsav);
        for (int i = 0; i < imax; ++i) {
            Real dsav = dlo + i * dstp;
            d[i] = std::pow(10.0e0_rt, dsav);
        }
    }

    // it does not work on all machines (for GPUs) broadcast to other
    // procs from managed memory.  So instead we'll read into a local
    // buffer, broadcast that, and then copy that into the managed
    // memory.

    amrex::Vector<Real> f_local(static_cast<size_t>(9) * imax * jmax);
    amrex::Vector<Real> dpdf_local(static_cast<size_t>(4) * imax * jmax);
    amrex::Vector<Real> ef_local(static_cast<size_t>(4) * imax * jmax);
    amrex::Vector<Real> xf_local(static_cast<size_t>(4) * imax * jmax);

    if (amrex::ParallelDescriptor::IOProcessor()) {

        // open the table
        std::ifstream table;
        table.open("helm_table.dat");

        if (!table.is_open()) {
            // the table was not present or we could not open it; abort
            amrex::Error("helm_table.dat could not be opened");
        }

        std::string line;

        // read in the free energy table
        int idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading free energy from helm_table.dat");
                }
                std::istringstream data(line);
                data >> f_local[idx] >> f_local[idx+3] >> f_local[idx+1]
                     >> f_local[idx+4] >> f_local[idx+2] >> f_local[idx+5]
                     >> f_local[idx+6] >> f_local[idx+7] >> f_local[idx+8];
                idx += 9;
            }
        }

        // read the pressure derivative with density table
        idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading pressure derivative from helm_table.dat");
                }
                std::istringstream data(line);
                data >> dpdf_local[idx] >> dpdf_local[idx+2]
                     >> dpdf_local[idx+1] >> dpdf_local[idx+3];
                idx += 4;
            }
        }

        // read the electron chemical potential table
        idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading electron chemical potential from helm_table.dat");
                }
                std::istringstream data(line);
                data >> ef_local[idx] >> ef_local[idx+2]
                     >> ef_local[idx+1] >> ef_local[idx+3];
                idx += 4;
            }
        }

        // read the number density table
        idx = 0;
        for (int j = 0; j < jmax; ++j) {
            for (int i = 0; i < imax; ++i) {
                std::getline(table, line);
                if (line.empty()) {
                    amrex::Error("Error reading number density from helm_table.dat");
                }
                std::istringstream data(line);
                data >> xf_local[idx] >> xf_local[idx+2]
                     >> xf_local[idx+1] >> xf_local[idx+3];
                idx += 4;
            }
        }

        table.close();

    }

    amrex::ParallelDescriptor::Bcast(f_local.data(),    static_cast<size_t>(9) * imax * jmax);
    amrex::ParallelDescriptor::Bcast(dpdf_local.data(), static_cast<size_t>(4) * imax * jmax);
    amrex::ParallelDescriptor::Bcast(ef_local.data(),   static_cast<size_t>(4) * imax * jmax);
    amrex::ParallelDescriptor::Bcast(xf_local.data(),   static_cast<size_t>(4) * imax * jmax);

    // now copy into managed memory
    int idx = 0;
    for (int j = 0; j < jmax; ++j) {  // NOLINT(modernize-loop-convert)
        for (int i = 0; i < imax; ++i) {
            for (int m = 0; m < 9; ++m) {
                f[j][i][m] = f_local[idx];
                idx++;
            }
        }
    }

    idx = 0;
    for (int j = 0; j < jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            for (int m = 0; m < 4; ++m) {
                dpdf[j][i][m] = dpdf_local[idx];
                ef[j][i][m] = ef_local[idx];
                xf[j][i][m] = xf_local[idx];
                idx++;
            }
        }
    }

    // construct the temperature and density deltas and their inverses
    for (int j = 0; j < jmax-1; ++j)
    {
        dth         = t[j+1] - t[j];
        dt2         = dth * dth;
        dti         = 1.0e0_rt / dth;
        dt2i        = 1.0e0_rt / dt2;
        dt_sav[j]   = dth;
        dt2_sav[j]  = dt2;
        dti_sav[j]  = dti;
        dt2i_sav[j] = dt2i;
    }

    for (int i = 0; i < imax-1; ++i)
    {
        dd          = d[i+1] - d[i];
        dd2         = dd * dd;
        ddi         = 1.0e0_rt / dd;
        dd2i        = 1.0e0_rt / dd2;
        dd_sav[i]   = dd;
        dd2_sav[i]  = dd2;
        ddi_sav[i]  = ddi;
        dd2i_sav[i] = dd2i;
    }

    // Set up the minimum and maximum possible densities.

    EOSData::mintemp = std::pow(10.e0_rt, tlo);
    EOSData::maxtemp = std::pow(10.e0_rt, thi);
    EOSData::mindens = std::pow(10.e0_rt, dlo);
    EOSData::maxdens = std::pow(10.e0_rt, dhi);
}



AMREX_INLINE
void actual_eos_finalize ()
{
}

#endif
