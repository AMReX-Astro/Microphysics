
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>

#include <nse_table.H>

void nse_cell_c()
{


    std::cout << "rho, T, Ye = " << density << " " << temperature << " " << ye << std::endl;

    // check the indices

    Real logrho = std::log10(density);
    Real logT = std::log10(temperature);

    int ir0 = nse_get_logrho_index(logrho);
    int it0 = nse_get_logT_index(logT);
    int ic0 = nse_get_ye_index(ye);

    std::cout << "density value brackets: "
              << nse_table_logrho(ir0) << " " << logrho << " " << nse_table_logrho(ir0+1) << std::endl;

    std::cout << "temperature value brackets: "
              << nse_table_logT(it0) << " " << logT << " " << nse_table_logT(it0+1) << std::endl;

    std::cout << "Ye value brackets: "
              << nse_table_ye(ic0) << " " << ye << " " << nse_table_ye(ic0+1) << std::endl;

    std::cout << std::endl;

    // now try 1-d interpolation in T -- for cubic interpolation, we
    // need to offset the indices one to the left

    ir0--;
    ir0 = amrex::max(1, amrex::min(table::nden-3, ir0));

    it0--;
    it0 = amrex::max(1, amrex::min(table::ntemp-3, it0));

    ic0--;
    ic0 = amrex::max(1, amrex::min(table::nye-3, ic0));


    Real rhos[] = {nse_table_logrho(ir0),
                   nse_table_logrho(ir0+1),
                   nse_table_logrho(ir0+2),
                   nse_table_logrho(ir0+3)};

    Real Ts[] = {nse_table_logT(it0),
                 nse_table_logT(it0+1),
                 nse_table_logT(it0+2),
                 nse_table_logT(it0+3)};

    Real yes[] = {nse_table_ye(ic0),
                  nse_table_ye(ic0+1),
                  nse_table_ye(ic0+2),
                  nse_table_ye(ic0+3)};


    {
        std::cout << "4 rho values: " << std::endl;

        Real _d[] = {table::abartab(nse_idx(ir0, it0, ic0)),
                     table::abartab(nse_idx(ir0+1, it0, ic0)),
                     table::abartab(nse_idx(ir0+2, it0, ic0)),
                     table::abartab(nse_idx(ir0+3, it0, ic0))};

        for (int ii = 0; ii < 4; ++ii) {
            std::cout << rhos[ii] << " " << _d[ii] << std::endl;
        }
        std::cout << std::endl;

        std::cout << "cubic interpolated value: " <<
            cubic(rhos, _d, table::dlogrho, logrho) << std::endl << std::endl;
    }

    {
        std::cout << "4 T values: " << std::endl;

        Real _d[] = {table::abartab(nse_idx(ir0, it0, ic0)),
                     table::abartab(nse_idx(ir0, it0+1, ic0)),
                     table::abartab(nse_idx(ir0, it0+2, ic0)),
                     table::abartab(nse_idx(ir0, it0+3, ic0))};

        for (int ii = 0; ii < 4; ++ii) {
            std::cout << Ts[ii] << " " << _d[ii] << std::endl;
        }
        std::cout << std::endl;

        std::cout << "cubic interpolated value: " <<
            cubic(Ts, _d, table::dlogT, logT) << std::endl << std::endl;

    }

    {
        std::cout << "4 Ye values: " << std::endl;

        Real _d[] = {table::abartab(nse_idx(ir0, it0, ic0)),
                     table::abartab(nse_idx(ir0, it0, ic0+1)),
                     table::abartab(nse_idx(ir0, it0, ic0+2)),
                     table::abartab(nse_idx(ir0, it0, ic0+3))};

        for (int ii = 0; ii < 4; ++ii) {
            std::cout << yes[ii] << " " << _d[ii] << std::endl;
        }
        std::cout << std::endl;

        std::cout << "cubic interpolated value: " <<
            cubic(yes, _d, -table::dye, ye) << std::endl << std::endl;

    }

    std::cout << "tricubic interpolated value: " <<
        tricubic(ir0, it0, ic0, logrho, logT, ye, table::abartab) << std::endl;

}
