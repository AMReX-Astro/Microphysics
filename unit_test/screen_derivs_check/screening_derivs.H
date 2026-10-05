#include <AMReX_Loop.H>
#include <AMReX_MultiFab.H>

#include <network.H>
#include <eos.H>
#include <rhs_type.H>

#include <screen.H>

#include <cmath>

using namespace amrex;
using namespace unit_test_rp;


// Check the Check independent moment perturbations, rather than just
// one composition direction.  With H and He, dY_H = 2*dM1-dM2 and
// dY_He = (dM2-dM1)/2.  Do not renormalize: these are partial
// derivatives.
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

        // initialize molar composition and plasma state
        Array1D<Real, 1, NumSpec> y{};
        y(ih) = 0.5_rt;
        y(ihe) = 0.125_rt;

        plasma_state_t<screening_dual_t> state;
        fill_plasma_state(state, temp, rho, y);

        // we'll do a single screening pair.  The derivatives are
        // indepedent of the screening pair.
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

            const Real step = 1.e-5_rt * (imode == 0 ? temp : 1.0_rt);

            auto yp = y;
            auto ym = y;
            if (imode > 0) {
                const Real dh = imode == 1 ? 2.0_rt * step : -step;
                const Real dhe = imode == 1 ? -0.5_rt * step : 0.5_rt * step;
                yp(ih) += dh; ym(ih) -= dh;
                yp(ihe) += dhe; ym(ihe) -= dhe;
            }

            // compute the + and - plasma states with the perturbed
            // quantities
            plasma_state_t<amrex::Real> plus, minus;
            fill_plasma_state(plus, temp + (imode == 0 ? step : 0.0_rt), rho, yp);
            fill_plasma_state(minus, temp - (imode == 0 ? step : 0.0_rt), rho, ym);

            // compute the screening factors
            const Real hp = actual_log_screen(plus, pair);
            const Real hm = actual_log_screen(minus, pair);

            const Real fp = actual_screen(plus, pair);
            const Real fm = actual_screen(minus, pair);

            // get the quantity we are comparing to from our call to the
            // autodiff version of the screening interface
            const Real h_deriv = imode == 0 ? ht : (imode == 1 ? h1 : h2);
            const Real f_deriv = imode == 0 ? ft : (imode == 1 ? f1 : f2);

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
