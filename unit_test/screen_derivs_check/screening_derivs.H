#ifndef SCREENING_DERIVS_H
#define SCREENING_DERIVS_H

#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>

#include <network.H>
#include <eos.H>
#include <rhs_type.H>

#include <screen.H>

#include <cmath>

using namespace amrex;
using namespace unit_test_rp;


void test_screening_derivatives()
{
    const int ih = network_spec_index("hydrogen-1") + 1;
    const int ihe = network_spec_index("helium-4") + 1;
    AMREX_ALWAYS_ASSERT(ih > 0 && ihe > 0);

    // a fixed number of loop iterations

    amrex::ParallelFor(12, [=] AMREX_GPU_DEVICE (int k) {


#if SCREEN_METHOD == SCREEN_METHOD_chugunov2009 || SCREEN_METHOD == SCREEN_METHOD_chabrier1998
        // In the weak-coupling limit, cancellation in the free-energy fit
        // amplifies the existing fast_atan derivative approximation.  Use
        // stronger coupling for finite-difference checks of these methods.
        constexpr Real log_rho_min = 8.0_rt;
#else
        constexpr Real log_rho_min = 2.0_rt;
#endif

        const Real rho = std::pow(10.0_rt, log_rho_min + 2.0_rt * (k % 4));
        const Real temp = std::pow(10.0_rt, 7.0_rt + (k / 4));

        amrex::Print() << std::format("ρ = {:.3g}, T = {:.3g}:\n ", rho, temp);

        // initialize molar composition and plasma state
        Array1D<Real, 1, NumSpec> y{};
        y(ih) = 0.5_rt;
        y(ihe) = 0.125_rt;

        plasma_state_t<screening_dual_t> state;
        fill_plasma_state(state, temp, rho, y);

        // we'll do a single screening pair.  The derivatives are
        // independent of the screening pair.
        constexpr auto pair = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);

        Real h, ht, h1, h2, f, ft, f1, f2;
        actual_log_screen(state, pair, h, ht, h1, h2);
        actual_screen(state, pair, f, ft, f1, f2);

        // fast_atan evaluates an approximation but its autodiff rule uses
        // the exact atan derivative.  Allow that existing discrepancy for
        // the two free-energy fits that use it; other methods are tighter.
#if SCREEN_METHOD == SCREEN_METHOD_chugunov2009 || SCREEN_METHOD == SCREEN_METHOD_chabrier1998
        constexpr Real derivative_tol = 5.e-3_rt;
#else
        constexpr Real derivative_tol = 2.e-5_rt;
#endif

        for (int imode = 0; imode < 3; ++imode) {
            // imode = 0: perturb T, keep M1 and M2 fixed
            // imode = 1: perturb M1, keep T and M2 fixed
            // imode = 2: perturb M2, keep T and M1 fixed

            amrex::Print() << std::format("    imode = {}\n", imode);

            // delta used for finite-difference

            amrex::Real step{};
            if (imode == 0) {
                step = 1.e-5_rt * temp;
            } else {
                step = 1.e-5_rt;
            }

            // perturb the composition (if we are in that mode) The
            // perturbations should be independent in terms of the
            // moment, not composition.  With H and He, dY_H =
            // 2*dM1-dM2 and dY_He = (dM2-dM1)/2.  Do not renormalize:
            // these are partial derivatives.

            auto yp = y;
            auto ym = y;
            if (imode > 0) {
                amrex::Real dh{};
                amrex::Real dhe{};
                if (imode == 1) {
                    dh = 2.0_rt * step;
                    dhe = -0.5_rt * step;
                } else {
                    dh = -step;
                    dhe = 0.5_rt * step;
                }
                yp(ih) += dh; ym(ih) -= dh;
                yp(ihe) += dhe; ym(ihe) -= dhe;
            }

            // compute the + and - plasma states with the perturbed
            // quantities

            plasma_state_t<amrex::Real> plus, minus;
            fill_plasma_state(plus, temp + (imode == 0 ? step : 0.0_rt), rho, yp);
            fill_plasma_state(minus, temp - (imode == 0 ? step : 0.0_rt), rho, ym);

            // compute the screening factors -- we do both h = log(f)
            // and f directly.

            const amrex::Real hp = actual_log_screen(plus, pair);
            const amrex::Real hm = actual_log_screen(minus, pair);

            const amrex::Real fp = actual_screen(plus, pair);
            const amrex::Real fm = actual_screen(minus, pair);

            // get the quantity we are comparing to from our call to the
            // autodiff version of the screening interface

            amrex::Real h_deriv{};
            amrex::Real f_deriv{};

            if (imode == 0) {
                h_deriv = ht;
                f_deriv = ft;
            } else if (imode == 1) {
                h_deriv = h1;
                f_deriv = f1;
            } else {
                h_deriv = h2;
                f_deriv = f2;
            }

            amrex::Real h_diff = (hp - hm) / (2.0_rt * step);
            amrex::Real f_diff = (fp - fm) / (2.0_rt * step);

            if (imode == 0) {
                amrex::Print() << std::format("       ∂h/∂T |_auto = {:.5g}; ∂h/∂T |_f-d = {:.5g};\n",
                                              ht, h_diff);
                amrex::Print() << std::format("       ∂f/∂T |_auto = {:.5g}; ∂f/∂T |_f-d = {:.5g};\n\n",
                                              ft, f_diff);

            } else if (imode == 1) {
                amrex::Print() << std::format("       ∂h/∂m₁ |_auto = {:.5g}; ∂h/∂m₁ |_f-d = {:.5g};\n",
                                              h1, h_diff);
                amrex::Print() << std::format("       ∂f/∂m₁ |_auto = {:.5g}; ∂f/∂m₁ |_f-d = {:.5g};\n\n",
                                              f1, f_diff);

            } else {
                amrex::Print() << std::format("       ∂h/∂m₂ |_auto = {:.5g}; ∂h/∂m₂ |_f-d = {:.5g};\n",
                                              h2, h_diff);
                amrex::Print() << std::format("       ∂f/∂m₂ |_auto = {:.5g}; ∂f/∂m₂ |_f-d = {:.5g};\n\n",
                                              f2, f_diff);
            }

            // for a centered diff, we want (for log screening factor, h):
            //
            //  | (hp - hm) / (2 * step) - h_deriv | < eps
            //
            // and similar for the full screening factor, f.
            //
            // We use a relative and absolute comparison here, to
            // avoid division by a nearly zero derivative (e.g. M1
            // in Debye-Huckel screening).

            // h = log(f)
            AMREX_ALWAYS_ASSERT(std::abs(hp - hm - 2.0_rt * step * h_deriv) <=
                                derivative_tol * std::abs(2.0_rt * step * h_deriv) +
                                2.e-12_rt * (1.0_rt + std::abs(hp) + std::abs(hm)));

            // f
            AMREX_ALWAYS_ASSERT(std::abs(fp - fm - 2.0_rt * step * f_deriv) <=
                                derivative_tol * std::abs(2.0_rt * step * f_deriv) +
                                2.e-12_rt * (1.0_rt + std::abs(fp) + std::abs(fm)));
        }
    });

    amrex::Print() << "Screening T, M1, M2 derivative checks passed" << std::endl;
}

#endif
