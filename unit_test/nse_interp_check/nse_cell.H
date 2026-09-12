#ifndef NSE_CELL_H
#define NSE_CELL_H

#include <cmath>
#include <format>

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <iostream>

#include <nse_table.H>
#include <nse_table_type.H>
#include <nse_eos.H>

AMREX_INLINE
void nse_cell_c()
{

    std::cout << std::endl;

    std::cout << "rho, T, Ye = "
              << unit_test_rp::density << " "
              << unit_test_rp::temperature << " "
              << unit_test_rp::ye << std::endl << std::endl;

    // check the indices

    amrex::Real logrho = std::log10(unit_test_rp::density);
    amrex::Real logT = std::log10(unit_test_rp::temperature);

    int ir0 = nse_get_logrho_index(logrho);
    int it0 = nse_get_logT_index(logT);
    int ic0 = nse_get_ye_index(unit_test_rp::ye);

    std::cout << "density value brackets: "
              << nse_table_logrho(ir0) << " < " << logrho << " < "
              << nse_table_logrho(ir0+1) << std::endl;

    std::cout << "temperature value brackets: "
              << nse_table_logT(it0) << " < " << logT << " < "
              << nse_table_logT(it0+1) << std::endl;

    std::cout << "Ye value brackets: "
              << nse_table_ye(ic0) << " < " << unit_test_rp::ye << " < "
              << nse_table_ye(ic0+1) << std::endl;

    std::cout << std::endl;

    // now try 1-d interpolation in rho, T, and Ye -- for cubic interpolation, we
    // need to offset the indices one to the left

    std::cout << "testing 1-d cubic interpolation at input state\n\n";

    {

        // density interpolation
        // note: assuming we are not at a boundary

        amrex::Real rhos[] = {nse_table_logrho(ir0-1),
                              nse_table_logrho(ir0),
                              nse_table_logrho(ir0+1),
                              nse_table_logrho(ir0+2)};

        std::cout << "    4 rho values: " << std::endl;

        amrex::Real _d[] = {nse_table::abartab(nse_idx(ir0-1, it0, ic0)),
                            nse_table::abartab(nse_idx(ir0, it0, ic0)),
                            nse_table::abartab(nse_idx(ir0+1, it0, ic0)),
                            nse_table::abartab(nse_idx(ir0+2, it0, ic0))};

        for (int ii = 0; ii < 4; ++ii) {
            std::cout << std::format("    {:5.3f} {:15.10g}\n", rhos[ii], _d[ii]);
        }
        std::cout << std::endl;

        std::cout << "    cubic interpolated value: " <<
            cubic(rhos, _d, nse_table_size::dlogrho, logrho, false, false) << std::endl << std::endl;
    }

    {

        // temperature interpolation

        amrex::Real Ts[] = {nse_table_logT(it0-1),
                            nse_table_logT(it0),
                            nse_table_logT(it0+1),
                            nse_table_logT(it0+2)};

        std::cout << "    4 T values: " << std::endl;

        amrex::Real _d[] = {nse_table::abartab(nse_idx(ir0, it0-1, ic0)),
                            nse_table::abartab(nse_idx(ir0, it0, ic0)),
                            nse_table::abartab(nse_idx(ir0, it0+1, ic0)),
                            nse_table::abartab(nse_idx(ir0, it0+2, ic0))};

        for (int ii = 0; ii < 4; ++ii) {
            std::cout << std::format("    {:5.3f} {:15.10g}\n", Ts[ii], _d[ii]);
        }
        std::cout << std::endl;

        std::cout << "    cubic interpolated value: " <<
            cubic(Ts, _d, nse_table_size::dlogT, logT, false, false) << std::endl << std::endl;

    }

    {

        // Ye interpolation

        amrex::Real yes[] = {nse_table_ye(ic0-1),
                             nse_table_ye(ic0),
                             nse_table_ye(ic0+1),
                             nse_table_ye(ic0+2)};

        std::cout << "    4 Ye values: " << std::endl;

        amrex::Real _d[] = {nse_table::abartab(nse_idx(ir0, it0, ic0-1)),
                            nse_table::abartab(nse_idx(ir0, it0, ic0)),
                            nse_table::abartab(nse_idx(ir0, it0, ic0+1)),
                            nse_table::abartab(nse_idx(ir0, it0, ic0+2))};

        for (int ii = 0; ii < 4; ++ii) {
            std::cout << std::format("    {:5.3f} {:15.10g}\n", yes[ii], _d[ii]);
        }
        std::cout << std::endl;

        std::cout << "    cubic interpolated value: " <<
            cubic(yes, _d, -nse_table_size::dye, unit_test_rp::ye, false, false) << std::endl << std::endl;

    }

    std::cout << std::endl;

    std::cout << "tricubic interpolated values at input state\n\n";

    nse_table_t nse_state;
    nse_state.T = unit_test_rp::temperature;
    nse_state.rho = unit_test_rp::density;
    nse_state.Ye = unit_test_rp::ye;

    nse_interp(nse_state);

    std::cout << "    abar    = " << nse_state.abar << std::endl;
    std::cout << "    bea     = " << nse_state.bea << std::endl;
    std::cout << "    dyedt   = " << nse_state.dyedt << std::endl;
    std::cout << "    dbeadt = " << nse_state.dbeadt << std::endl;
    std::cout << "    e_nu    = " << nse_state.e_nu << std::endl;
    std::cout << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "    X(" << short_spec_names_cxx[n] << ") = " << nse_state.X[n] << std::endl;
    }


    // temp derivatives

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "testing temperature derivatives of cubic" << std::endl << std::endl;

    std::cout << "    first finite-difference derivatives\n";

    nse_state.T = unit_test_rp::temperature;
    nse_state.rho = unit_test_rp::density;
    nse_state.Ye = unit_test_rp::ye;

    nse_interp(nse_state);

    amrex::Real abar_old = nse_state.abar;
    amrex::Real bea_old = nse_state.bea;
    amrex::Real T_old = nse_state.T;

    const amrex::Real eps = 1.e-8_rt;
    nse_state.T *= (1.0_rt + eps);

    nse_interp(nse_state);

    std::cout << "        dAbar/dT = " << (nse_state.abar - abar_old) / (nse_state.T - T_old) << std::endl;
    std::cout << "        dbea/dT = " << (nse_state.bea - bea_old) / (nse_state.T - T_old) << std::endl;

    std::cout << std::endl;
    std::cout << "    now using derivative of the interpolant" << std::endl;

    amrex::Real dabardT = nse_interp_dT(unit_test_rp::temperature, unit_test_rp::density, unit_test_rp::ye,
                                 nse_table::abartab);

    amrex::Real dbeadT = nse_interp_dT(unit_test_rp::temperature, unit_test_rp::density, unit_test_rp::ye,
                                 nse_table::beatab);

    std::cout << "        dAbar/dT = " << dabardT << std::endl;
    std::cout << "        dbea/dT = " << dbeadT << std::endl;

    std::cout << std::endl;

    // dens derivatives

    std::cout << "testing density derivatives of cubic" << std::endl << std::endl;

    std::cout << "    first finite-difference derivatives" << std::endl;

    nse_state.T = unit_test_rp::temperature;
    nse_state.rho = unit_test_rp::density;
    nse_state.Ye = unit_test_rp::ye;

    nse_interp(nse_state);

    amrex::Real rho_old = nse_state.rho;

    nse_state.rho *= (1.0_rt + eps);

    nse_interp(nse_state);

    std::cout << "        dAbar/drho = " << (nse_state.abar - abar_old) / (nse_state.rho - rho_old) << std::endl;
    std::cout << "        dbea/drho = " << (nse_state.bea - bea_old) / (nse_state.rho - rho_old) << std::endl;

    std::cout << std::endl;
    std::cout << "    now using derivative of the interpolant" << std::endl;

    amrex::Real dabardrho = nse_interp_drho(unit_test_rp::temperature, unit_test_rp::density, unit_test_rp::ye,
                                     nse_table::abartab);

    amrex::Real dbeadrho = nse_interp_drho(unit_test_rp::temperature, unit_test_rp::density, unit_test_rp::ye,
                                    nse_table::beatab);

    std::cout << "        dAbar/drho = " << dabardrho << std::endl;
    std::cout << "        dbea/drho = " << dbeadrho << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;

    //
    // EOS testing
    //

    // now we test the EOS, in particular, we want to ensure that
    // given e, rho, Y_e, we can find a T that is consistent with both
    // the EOS and the NSE table

    // attempt 1: using the normal EOS interfaces
    //
    // we'll do:
    //
    // T, rho, Ye -> Abar
    // T, rho, Ye, Abar -> e
    //
    // then perturb e to e' and do:
    //
    // e', rho, Ye -> T'
    // T', rho, Ye -> Abar'
    //
    // and finally check
    //
    // T', rho, Ye, Abar' -> e' ???

    {

        std::cout << "EOS T from e consistency check (e should remain unchanged)" << std::endl << std::endl;

        // first get the abar consistent with our inputs and find e

        nse_state.T = unit_test_rp::temperature;
        nse_state.rho = unit_test_rp::density;
        nse_state.Ye = unit_test_rp::ye;

        nse_interp(nse_state);

        amrex::Real abar_orig = nse_state.abar;

        eos_t eos_state;

        eos_state.T = unit_test_rp::temperature;
        eos_state.rho = unit_test_rp::density;
        eos_state.aux[iye] = nse_state.Ye;
        eos_state.aux[iabar] = nse_state.abar;

        eos(eos_input_rt, eos_state);

        amrex::Real e_orig = eos_state.e;

        // now perturb e and find T, then redo NSE to get abar and finally
        // see if we get back our new e

        amrex::Real e_new = eos_state.e * 1.05;

        eos_state.e = e_new;
        eos(eos_input_re, eos_state);

        nse_state.T = eos_state.T;
        nse_interp(nse_state);

        eos_state.aux[iabar] = nse_state.abar;
        eos(eos_input_rt, eos_state);

        std::cout << "    old method: invert EOS without considering NSE changes:" << std::endl;
        std::cout << "        change in e: " << eos_state.e << " " << e_new << std::endl;

        // attempt 2:
        // now we try the new interface.  This effectively does:
        // e', rho, Ye -> Abar', T'

        eos_state.T = unit_test_rp::temperature;
        eos_state.e = e_new;
        eos_state.rho = unit_test_rp::density;
        eos_state.aux[iye] = unit_test_rp::ye;
        eos_state.aux[iabar] = abar_orig;

        amrex::Real abar_start = eos_state.aux[iabar];

        nse_T_abar_from_e(eos_state.rho, eos_state.e, eos_state.aux[iye],
                          eos_state.T, eos_state.aux[iabar]);

        std::cout << std::endl;
        std::cout << "    NSE method: use the nse_T_abar_from_e() function:" << std::endl;
        std::cout << "        updated T: " << eos_state.T << std::endl;
        std::cout << "        change in abar: " << abar_start << " " << eos_state.aux[iabar] << std::endl;

        // now check if we get back the correct e!

        eos(eos_input_rt, eos_state);

        std::cout << "        change in e: " << eos_state.e << " " << e_new << std::endl;
        std::cout << std::endl;
    }


    {

        // now redo it for pressure

        std::cout << "EOS T from p consistency check (p should remain unchanged)" << std::endl << std::endl;

        nse_state.T = unit_test_rp::temperature;
        nse_state.rho = unit_test_rp::density;
        nse_state.Ye = unit_test_rp::ye;

        nse_interp(nse_state);

        amrex::Real abar_orig = nse_state.abar;

        eos_t eos_state;

        eos_state.T = unit_test_rp::temperature;
        eos_state.rho = unit_test_rp::density;
        eos_state.aux[iye] = nse_state.Ye;
        eos_state.aux[iabar] = nse_state.abar;

        eos(eos_input_rt, eos_state);

        amrex::Real p_orig = eos_state.p;

        // now perturb p and find T, then redo NSE to get abar and finally
        // see if we get back our new p

        amrex::Real p_new = eos_state.p * 1.05;

        eos_state.p = p_new;
        eos(eos_input_rp, eos_state);

        nse_state.T = eos_state.T;
        nse_interp(nse_state);

        eos_state.aux[iabar] = nse_state.abar;
        eos(eos_input_rt, eos_state);

        std::cout << "    old method: invert EOS without considering NSE changes:" << std::endl;
        std::cout << "        change in p: " << eos_state.p << " " << p_new << std::endl;

        // attempt 2:
        // now we try the new interface.  This effectively does:
        // p', rho, Ye -> Abar', T'

        eos_state.T = unit_test_rp::temperature;
        eos_state.p = p_new;
        eos_state.rho = unit_test_rp::density;
        eos_state.aux[iye] = unit_test_rp::ye;
        eos_state.aux[iabar] = abar_orig;

        amrex::Real abar_start = eos_state.aux[iabar];

        nse_T_abar_from_p(eos_state.rho, eos_state.p, eos_state.aux[iye],
                          eos_state.T, eos_state.aux[iabar]);

        std::cout << std::endl;
        std::cout << "    NSE method: use the nse_T_abar_from_p() function:" << std::endl;
        std::cout << "        updated T: " << eos_state.T << std::endl;
        std::cout << "        change in abar: " << abar_start << " " << eos_state.aux[iabar] << std::endl;

        // now check if we get back the correct p!

        eos(eos_input_rt, eos_state);

        std::cout << "        change in p: " << eos_state.p << " " << p_new << std::endl;
        std::cout << std::endl;
    }


    // now we test the EOS inversion for finding rho given e or p.  The idea
    // is the same as above.  We don't have an eos_input_te, so we will only
    // check p

    {

        std::cout << "EOS rho from p consistency check (p should remain unchanged)" << std::endl << std::endl;

        nse_state.T = unit_test_rp::temperature;
        nse_state.rho = unit_test_rp::density;
        nse_state.Ye = unit_test_rp::ye;

        nse_interp(nse_state);

        amrex::Real abar_orig = nse_state.abar;

        eos_t eos_state;

        eos_state.T = unit_test_rp::temperature;
        eos_state.rho = unit_test_rp::density;
        eos_state.aux[iye] = nse_state.Ye;
        eos_state.aux[iabar] = nse_state.abar;

        eos(eos_input_rt, eos_state);

        amrex::Real p_orig = eos_state.p;

        // now perturb p and find rho, then redo NSE to get abar and finally
        // see if we get back our new p

        amrex::Real p_new = eos_state.p * 1.05;

        eos_state.p = p_new;
        eos(eos_input_tp, eos_state);

        nse_state.rho = eos_state.rho;
        nse_interp(nse_state);

        eos_state.aux[iabar] = nse_state.abar;
        eos(eos_input_rt, eos_state);

        std::cout << "    old method: invert EOS without considering NSE changes:" << std::endl;
        std::cout << "        change in p: " << eos_state.p << " " << p_new << std::endl;

        // attempt 2:
        // now we try the new interface.  This effectively does:
        // p', T, Ye -> Abar', rho'

        eos_state.T = unit_test_rp::temperature;
        eos_state.p = p_new;
        eos_state.rho = unit_test_rp::density;
        eos_state.aux[iye] = unit_test_rp::ye;
        eos_state.aux[iabar] = abar_orig;

        amrex::Real abar_start = eos_state.aux[iabar];

        nse_rho_abar_from_p(eos_state.T, eos_state.p, eos_state.aux[iye],
                            eos_state.rho, eos_state.aux[iabar]);

        std::cout << std::endl;
        std::cout << "    NSE method: use the nse_rho_abar_from_p() function:" << std::endl;

        std::cout << "        updated T: " << eos_state.T << std::endl;
        std::cout << "        change in abar: " << abar_start << " " << eos_state.aux[iabar] << std::endl;

        // now check if we get back the correct p!

        eos(eos_input_rt, eos_state);

        std::cout << "        change in p: " << eos_state.p << " " << p_new << std::endl;
        std::cout << std::endl;
    }

    {

        // now explore a troublesome state where tricubic
        // interpolation can lead to changing the sign of the neutrino
        // energy.

        const amrex::Real T = 3667943011.533362;
        const amrex::Real rho = 9434630294.888599;
        const amrex::Real Ye = 0.433159096449449;

        amrex::Real logrho = std::log10(rho);
        amrex::Real logT = std::log10(T);

        // shift the indices such that it is the starting value
        // for the 4x4x4 cube we use for tricubic
        int ir0 = nse_get_logrho_index(logrho) - 1;
        ir0 = amrex::Clamp(ir0, 1, nse_table_size::nden-3);

        int it0 = nse_get_logT_index(logT) - 1;
        it0 = amrex::Clamp(it0, 1, nse_table_size::ntemp-3);

        int ic0 = nse_get_ye_index(Ye) - 1;
        ic0 = amrex::Clamp(ic0, 1, nse_table_size::nye-3);

        std::cout << std::endl;
        std::cout << "explore a troublesome state" << std::endl << std::endl;

        std::cout << std::format("    rho = {:20.15e}, T = {:20.15e}, Ye = {}\n",
                                 rho, T, Ye);
        std::cout << std::format("    indices of lowest corner of 4x4x4 cube: ir0 = {}, it0 = {}, ic0 = {}\n",
                                 ir0, it0, ic0);

        nse_table_t nse_state;
        nse_state.T = T;
        nse_state.rho = rho;
        nse_state.Ye = Ye;

        nse_interp(nse_state);

        std::cout << std::endl;
        std::cout << "    cubic interpolation gives: e_nu = " << nse_state.e_nu << std::endl;

        amrex::Real enu_min = std::numeric_limits<amrex::Real>::max();
        amrex::Real enu_max = std::numeric_limits<amrex::Real>::lowest();

        for (int ir = ir0; ir < ir0 + 4; ++ir) {
            for (int jt = it0; jt < it0 + 4; ++jt) {
                for (int ky = ic0; ky < ic0 + 4; ++ky) {

                    auto _enu = nse_table::enutab(nse_idx(ir, jt, ky));

                    enu_min = std::min(enu_min, _enu);
                    enu_max = std::max(enu_max, _enu);

                }
            }
        }

        std::cout << "    bounds of e_nu in the 64 cells used for tricubic: " << enu_min << " " << enu_max << std::endl;


        std::cout << std::endl;

    }




}
#endif
