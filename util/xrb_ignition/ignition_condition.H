#ifndef IGNITION_CONDITION_H
#define IGNITION_CONDITION_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <cassert>
#include <actual_conductivity.H>
#include <screen.H>
#include <microphysics_autodiff.H>
#include <unit_test.H>

/* mass fractions of CNO elements (metal) based on beta-decay timescale
   O14 -> N14 (70.641 s) and O15 -> N15 (122.24 s).*/

constexpr amrex::Real O14 = 0.352;
constexpr amrex::Real O15 = 0.648;

struct State
{
    // Radius and gravitational acceleration
    amrex::Real R, g;

    // Initial abundance for hydrogen, helium, and metal
    amrex::Real X0, Y0, Z0;

    // Local Accretion rate
    amrex::Real mdot;

    // Depletion Column Depth
    amrex::Real yd;

    // Top boundary column depth, temperature and flux
    // This is used as the initial condition from integration.
    amrex::Real yt, Tt, Ft;

    // Bot boundary column depth, temperature and flux
    // This is assumed to be the ignition condition
    amrex::Real yb, Tb, Fb;

    // Height below the top layer corresponding to yt,
    amrex::Real z;
};


template<int NEQ>
AMREX_FORCE_INLINE
void rk4_step(amrex::Real t, amrex::Real dt,
              const amrex::Real (&f)[NEQ], const State& state,
              void (*rhs)(amrex::Real, const amrex::Real (&)[NEQ],
                          const State&, amrex::Real (&)[NEQ]),
              amrex::Real (&fout)[NEQ])
{
    amrex::Real k1[NEQ], k2[NEQ], k3[NEQ], k4[NEQ];
    amrex::Real ftmp[NEQ];

    rhs(t, f, state, k1);

    for (int n = 0; n < NEQ; ++n)
        ftmp[n] = f[n] + 0.5_rt * dt * k1[n];
    rhs(t + 0.5_rt*dt, ftmp, state, k2);

    for (int n = 0; n < NEQ; ++n)
        ftmp[n] = f[n] + 0.5_rt * dt * k2[n];
    rhs(t + 0.5_rt*dt, ftmp, state, k3);

    for (int n = 0; n < NEQ; ++n)
        ftmp[n] = f[n] + dt * k3[n];
    rhs(t + dt, ftmp, state, k4);

    for (int n = 0; n < NEQ; ++n)
        fout[n] = f[n] + (dt / 6.0_rt) *
            (k1[n] + 2.0_rt*k2[n] + 2.0_rt*k3[n] + k4[n]);
}


template<int NEQ>
void rk4_adaptive(amrex::Real t0, amrex::Real tmax,
                  amrex::Real dt_init, amrex::Real eps,
                  const amrex::Real (&f0)[NEQ], const State& state,
                  void (*rhs)(amrex::Real, const amrex::Real (&)[NEQ],
                              const State&, amrex::Real (&)[NEQ]),
                  amrex::Real (&f_new)[NEQ])
{
    constexpr amrex::Real SAFETY = 0.9_rt;
    constexpr amrex::Real GROW   = -0.2_rt;   // 1/(order+1)
    constexpr amrex::Real SHRINK = -0.25_rt;  // 1/order
    constexpr amrex::Real TINY   = 1.e-30_rt;

    amrex::Real t = t0;
    amrex::Real dt = dt_init;

    for (int n = 0; n < NEQ; ++n) {
        f_new[n] = f0[n];
    }

    amrex::Real y_big[NEQ], y_half[NEQ], y_tmp[NEQ], yscal[NEQ];

    while (t < tmax) {

        if (t + dt > tmax) {
            dt = tmax - t;
        }

        // scaling for error norm
        for (int n = 0; n < NEQ; ++n) {
            yscal[n] = std::abs(f_new[n]) + TINY;
        }

        amrex::Real errmax;

        for (;;) {
            // one big step
            rk4_step<NEQ>(t, dt, f_new, state, rhs, y_big);

            // two half steps
            rk4_step<NEQ>(t, 0.5_rt*dt, f_new, state, rhs, y_tmp);
            rk4_step<NEQ>(t + 0.5_rt*dt, 0.5_rt*dt, y_tmp, state, rhs, y_half);

            // compute error
            errmax = 0.0_rt;
            for (int n = 0; n < NEQ; ++n) {
                amrex::Real err = std::abs(y_half[n] - y_big[n]) / yscal[n];
                errmax = std::max(errmax, err);
            }
            errmax /= eps;

            if (errmax <= 1.0_rt) {
                break;
            }

            // reduce timestep
            amrex::Real dt_new = SAFETY * dt * std::pow(errmax, SHRINK);
            dt = std::max(dt_new, 0.1_rt * dt);

            if (dt <= TINY) {
                amrex::Abort("rk4_adaptive: timestep underflow");
            }
        }

        // accept step
        t += dt;
        for (int n = 0; n < NEQ; ++n)
            f_new[n] = y_half[n];

        // increase timestep
        amrex::Real dt_new;
        if (errmax > 1.e-4_rt) {
            dt_new = SAFETY * dt * std::pow(errmax, GROW);
        } else {
            dt_new = 5.0_rt * dt;
        }
        dt = std::min(dt_new, tmax - t);
    }
}



template<int NEQ>
void rk4(const amrex::Real t0, const amrex::Real tmax, const amrex::Real dt_init,
         const amrex::Real (&f0)[NEQ], const State& state,
         void (*rhs)(const amrex::Real, const amrex::Real (&)[NEQ],
                     const State&, amrex::Real (&)[NEQ]),
         amrex::Real (&f_new)[NEQ]) {
    /* Runge-Kutta 4 IVP solver for multiple variables: df1/dt, df2dt, ...
       where NEQ is the number of variables. It just solves to tmax.

       Parameters
       ----------
       t0: Initial time
       tmax: Final time
       dt_init: initial Stepsize used for integration
       f0: Initial value corresponding to t0, f(t0) = f0.
       rhs: Function that computes the derivative, df / dt
     */

    // Set initial condition

    amrex::Real dt = dt_init;
    amrex::Real t = t0;

    for (int n = 0; n < NEQ; ++n) {
        f_new[n] = f0[n];
    }

    amrex::Real k1[NEQ], k2[NEQ], k3[NEQ], k4[NEQ];
    amrex::Real ftmp[NEQ], f[NEQ];

    while (t < tmax) {

        // Change stepsize so that we don't overshoot
        if (t + dt > tmax) {
            dt = tmax - t;
        }

        // Set current f
        for (int n = 0; n < NEQ; ++n) {
            f[n] = f_new[n];
        }

        rhs(t, f, state, k1);

        for (int n = 0; n < NEQ; ++n) {
            ftmp[n] = f[n] + 0.5_rt * dt * k1[n];
        }
        rhs(t + 0.5_rt * dt, ftmp, state, k2);

        for (int n = 0; n < NEQ; ++n) {
            ftmp[n] = f[n] + 0.5_rt * dt * k2[n];
        }
        rhs(t + 0.5_rt * dt, ftmp, state, k3);

        for (int n = 0; n < NEQ; ++n) {
            ftmp[n] = f[n] + dt * k3[n];
        }
        rhs(t + dt, ftmp, state, k4);

        for (int n = 0; n < NEQ; ++n) {
            f_new[n] = f[n] + (dt / 6.0_rt) *
                (k1[n] + 2.0_rt*k2[n] + 2.0_rt*k3[n] + k4[n]);
        }

        // update time
        t += dt;
    }
}


amrex::Real zbrent(amrex::Real a,
                   amrex::Real b,
                   amrex::Real tol,
                   int max_iter,
                   amrex::Real (*fcn)(amrex::Real, State&),
                   State& state) {
    /*
      Brent's method for finding root without Jacobian within interval [a, b]
     */

    // Compute the constraint equation at the boundary
    amrex::Real fa = fcn(a, state);
    amrex::Real fb = fcn(b, state);

    // Makre sure roots are bracketed
    // i.e. they have opposite sign
    if (fa * fb > 0.0_rt) {
        std::cout << "Warning: root not bracketed" << std::endl;
    }

    // Make large value on the left?
    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    amrex::Real c = a;
    amrex::Real fc = fa;

    // Flag to see if we went through bisection method previously
    bool mflag = true;

    // s is the new solution
    amrex::Real s = b;
    amrex::Real d = 0.0_rt;

    for (int iter = 0; iter < max_iter; ++iter) {

        if (fa != fc && fb != fc) {
            // Inverse quadratic interpolation if we have 3 distinct points
            // a=x_n-2, b=x_n-1, c=x_n
            s = a * fb * fc / ((fa - fb) * (fa - fc)) +
                b * fa * fc / ((fb - fa) * (fb - fc)) +
                c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            // Secant method otherwise
            s = b - fb * (b - a) / (fb - fa);
        }

        // Conditions to fall back to bisection
        if ((s < (3*a + b) * 0.25_rt || s > b) ||
            (mflag && std::abs(s - b) >= 0.5_rt * std::abs(b - c)) ||
            (!mflag && std::abs(s - b) >= 0.5_rt * std::abs(c - d)) ||
            (mflag && std::abs(b - c) < tol) ||
            (!mflag && std::abs(c - d) < tol)) {

            s = 0.5_rt * (a + b);
            mflag = true;
        } else {
            mflag = false;
        }

        amrex::Real fs = fcn(s, state);

        // d is previous value of c
        d = c;

        // c is previous value of b
        c = b;
        fc = fb;

        if (fa * fs < 0.0_rt) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }

        // Convergence test
        if (std::abs(b - a) < tol || std::abs(fb) < tol) {
            return b;
        }
    }

    amrex::Abort("Brent: did not converge");
    return b; // unreachable
}


void derivs(const amrex::Real y, const amrex::Real (&f)[3],
            const State& state, amrex::Real (&dfdy)[3]) {
    /*
      derivatives for integration:
      dfdy = (dT/dy, dF/dy, dz/dy) and f = (T, F, z)

      dT / dy = 3 κ / (4 a c) F / T^3
      dT / dy = F / (ϱ k_th)
      where k_th = (4 a c T^3) / (3 κ ϱ)

      dF/dy = -epsilon_H if X > 0, X = X0(1 - y/yd)

      dz/dy = -1 / ϱ
     */
    eos_t eos_state;
    eos_state.T = std::max(f[0], 1.0e7);
    eos_state.p = state.g * y;

    amrex::Real X = std::max(state.X0 * (1.0_rt - y / state.yd), 0.0_rt);
    eos_state.xn[Species::H1-1] = X;
    eos_state.xn[Species::He4-1] = 1.0_rt - X - state.Z0;
    eos_state.xn[Species::O14-1] = O14 * state.Z0;
    eos_state.xn[Species::O15-1] = O15 * state.Z0;

    // Find density via Temperature and pressure
    eos(eos_input_tp, eos_state);

    // Find conductivity
    actual_conductivity(eos_state);

    amrex::Real dTdy = f[1] / (eos_state.rho * eos_state.conductivity);

    dfdy[0] = dTdy;

    // Hot CNO burning
    amrex::Real dFdy = 0.0_rt;
    amrex::Real epsilon_H = 5.8e13 * (state.Z0 / 0.01);
    if (X > 0.0_rt) {
        dFdy = -epsilon_H;
    }
    dfdy[1] = dFdy;

    // Hydrostatic balance
    dfdy[2] = - 1.0_rt / eos_state.rho;
}


amrex::Real heating_term(const amrex::Real T0, const eos_t& eos_state) {
    /*
      Evaluates the heating term -- triple alpha process
      Use general fitting:
        ϵ_3⍺ = 5.3e21 F_scn rho5^2 Y^3 / T8^3 exp(-44 / T8)

      Also apply optional enhancement factor due to quick
      subsequent reaction rate: C12(p,g)N13(p,g)O14, where

        F_enh = 1 + Q12 / Q3a ≅ 1.9

      only when X(H) > 1.e-3.

      Note that we technically want temperature derivative.
      Although we technically have analytic form, but since
      the EOS is e(T,p), then ϱ is also temperature dependent.
      So to capture the effect, let's do a numerical difference.
      Use logarithmnic diff instead of linear diff for better accuracy.
     */

    // Q-value of C12 + 2p -> O14 ( C12(p,g)N13(p,g)O14 ) = 6.57 MeV
    // Q-value 3 He4 -> C12 = 7.274 MeV
    constexpr amrex::Real Q12 = 6.57_rt;
    constexpr amrex::Real Q3a = 7.274_rt;
    constexpr amrex::Real enhance_fac = 1.0 + Q12 / Q3a;

    // Call EOS to find density.
    auto tmp_state = eos_state;
    tmp_state.T = T0;
    eos(eos_input_tp, tmp_state);

    // Find screening factor

    constexpr int do_T_derivatives = 1;
    using number_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    number_t temp = tmp_state.T;
    if constexpr(do_T_derivatives) {
        autodiff::seed(temp);
    }

    plasma_state_t<number_t> pstate{};
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int n = 0; n < NumSpec; ++n) {
        Y(n+1) = tmp_state.xn[n] * aion_inv[n];
    }

    amrex::Real scor, scor2, dscor_dt, dscor2_dt;

    fill_plasma_state(pstate, temp, tmp_state.rho, Y);

    constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
    actual_screen(pstate, scn_fac, scor, dscor_dt);

    constexpr auto scn_fac2 = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
    actual_screen(pstate, scn_fac2, scor2, dscor2_dt);

    amrex::Real Fscn = scor * scor2;
    // amrex::Real dFscn_dT = scor * dscor2_dt + dscor_dt * scor2;

    // eps_3a = 5.3e21 Fscn rho5^2 Y^3 / T8^3 exp(-44 / T8)
    amrex::Real heat = 5.3e11 * Fscn * std::pow(tmp_state.rho, 2) *
        std::pow(tmp_state.xn[Species::He4-1], 3) * std::pow(tmp_state.T * 1.0e-8, -3) *
        std::exp(-44.0 * 1.0e8 / tmp_state.T);

    // If there is enough hydrogen for the rapid proton capture reactions:
    // C12(p,g)N13(p,g)O14, i.e. when X(H1) > 1.0e-3,
    // then apply an additional enhancement factor due to enhancement, ~ 1.9
    if (tmp_state.xn[Species::H1-1] > 1.0e-3_rt) {
        heat *= enhance_fac;
    }

    // Heating term is temperature derivative of epsilon_3a including screening enhancement.
    // Also note that since the EOS is e(T, p), this means ϱ is also temperature dependent.
    // Need dϱ/dT, but eos doesn't seem to provide this.
    // Let's use numerical jacobian instead of analytic, like settle version.

    // amrex::Real heat = (-3.0_rt * Fscn / tmp_state.T +
    //                     44.0_rt * Fscn / (tmp_state.T * tmp_state.T * 1.0e-8_rt) +
    //                     dFscn_dT) * epsilon_3a;

    return heat;
}


amrex::Real cooling_term(const amrex::Real T0, const amrex::Real yb,
                         const eos_t& eos_state) {
    /*
      Evaluates the cooling term.
      Get cooling rate by combining:

      dT/dy = 3 κ / (4 a c) F / T^3 and ε = -dF/dy

      ε ≅ a c / (3 κ) 1 / (y^2) T^4

      with k_th = (4 a c T^3) / (3 κ ϱ)

      ε ≅ 0.25 k_th ϱ T / y^2
      Note that since opacity is also temperature dependent,
      we should do a finite difference here.
     */

    // Call EOS to find conductivity
    auto tmp_state = eos_state;
    tmp_state.T = T0;
    eos(eos_input_tp, tmp_state);

    actual_conductivity(tmp_state);
    amrex::Real k_th = tmp_state.conductivity;
    amrex::Real cool = 0.25_rt * tmp_state.rho * tmp_state.T * k_th / (yb * yb);

    return cool;
}


amrex::Real yb_constraint_eq(amrex::Real log10yb_0, State& state) {
    /* Constraint equation to find ignition column depth
       given the guess ignition column depth, yb_0.
       Note the guess value is given in log10(yb) */

    state.yb = std::pow(10, log10yb_0);

    // Set upper boundary conditions, i.e. initial conditions: yt, Tt, and Ft

    // Pick arbitrary top boundary column depth
    // Follow "settle" version on Github
    state.yt = 1.e3_rt;

    // Compute the top boundary flux
    // Ft = F_b + epsilon_H * (yb or yd)
    // whichever is smaller.
    // If yb < yd, the ignition happens during hydrogen-helium layer
    // If yb > yd, the ignition happens during pure helium layer

    amrex::Real epsilon_H = 5.8e13 * (state.Z0 / 0.01);
    state.Ft = state.Fb + epsilon_H * std::min(state.yb, state.yd);

    // Compute the top boundary temperature
    // Given by the analytic radiative zero solution for constant flux
    // Tt = std::pow(3 κ Ft yt/ [a c], 1/4)

    // Here 2650 is 3 κ / (a c), which is precomputed.
    // where they assumed the κ at top is 0.22

    // If we were to compute kappa manually, we need eos.
    // i.e. we need two inputs. We have p=g * yt, so we need rho?

    state.Tt = std::pow(2650.0 * state.Ft * state.yt, 0.25);

    // Now set the lower boundary conditions: yb, Tb, and Fb
    // Here Fb is assumed to be given. Note that the original version
    // includes a special term for compressional heating, if thats included
    // then we have to solve for Fb. However, it is sufficient to include
    // an extra constant heating term to account for compressional heating.
    // yb is what we're after, so we need Tb.

    // The original version seems to use an ODE solver for this.
    // The idea is that opacity is also temperature dependent,
    // so we can't integrate dT/dy = 3 κ / (4 a c) F / T^3 easily
    // Also note that dT/dy depends on F, so we need to evolve F as well with dF/dy
    // dF/dy = -epsilon_H, where
    // Use RK4 to integrate from Tt and Ft all the way down to find T(yb) and F(yb)
    // How do I make sure F(yb) = Fb?
    // Note that we used Fb to find Ft: Ft = epsilon_H

    constexpr int NEQ = 3;
    amrex::Real f_array[NEQ];
    amrex::Real f0_array[NEQ] = {state.Tt, state.Ft, 0.0_rt};

    // t0=yt, tmax=yb, step_size=yt, f0={Tt, Ft, 0.0}
    // Integrates to find f_array corresponding to yb.
    rk4_adaptive<NEQ>(state.yt, state.yb, state.yt, 1.0e-6,
                      f0_array, state, derivs,
                      f_array);

    // Update solution Tb to State
    state.Tb = f_array[0];

    // Make sure solution to F(yb) is the same as Fb

    assert(std::abs(state.Fb - f_array[1]) < 1.0e-3 * state.Fb);

    // update relative height
    state.z = f_array[2];

    // Since we need density and opacity, get eos first.

    eos_t eos_state;
    eos_state.T = state.Tb;
    eos_state.p = state.g * state.yb;

    amrex::Real X = std::max(state.X0 * (1.0_rt - state.yb / state.yd), 0.0_rt);
    eos_state.xn[Species::H1-1] = X;
    eos_state.xn[Species::He4-1] = 1 - X - state.Z0;
    eos_state.xn[Species::O14-1] = O14 * state.Z0;
    eos_state.xn[Species::O15-1] = O15 * state.Z0;

    // After Tb is solved, assume we have runaway when
    // heating rate is larger than cooling rate, i.e.
    // d eps_3a / dT = d eps_cool / dT

    // To find temperature derivative of heating (triple-alpha) and cooling term
    // Do a logarithmic difference: de/dT = (e / T) dln(e)/dln(T)
    // dln(e)/dln(T) = [ln(e2) - ln(e1)] / [ln(T2) - ln(T1)] = ln(e2/e1) / ln(T2/T1)
    // This is more accurate than linear difference for expressions that involve exponentials.

    constexpr amrex::Real PERTURB_RATIO = 1.0e-3_rt;
    amrex::Real dT = state.Tb * PERTURB_RATIO;

    amrex::Real heat = heating_term(state.Tb, eos_state);
    amrex::Real heatp = heating_term(state.Tb + dT, eos_state);
    amrex::Real dheat_dT = heat * std::log(heatp/heat) /
        (state.Tb * std::log(1 + PERTURB_RATIO));

    amrex::Real cool = cooling_term(state.Tb, state.yb, eos_state);
    amrex::Real coolp = cooling_term(state.Tb + dT, state.yb, eos_state);
    amrex::Real dcool_dT = cool * std::log(coolp/cool) /
        (state.Tb * std::log(1 + PERTURB_RATIO));

    amrex::Real diff = (dheat_dT - dcool_dT) / dcool_dT;

    std::cout << std::scientific;
    std::cout << "heat rate is = " << dheat_dT << std::endl;
    std::cout << "cool rate is = " << dcool_dT << std::endl;
    std::cout << "yb = " << state.yb
              << ", Tb = " << state.Tb
              << ", (heat-cool) / cool = " << diff
              << std::endl;

    return diff;
}


void find_yb(State& state) {
    // Root-find using the constraint equation to find yb
    // Use zbrent method so that Jacobian is not required.
    // Assumes that yb is between 10^8 ~ 10^10

    zbrent(8.0_rt, 10.0_rt, 1.0e-3_rt, 1000, yb_constraint_eq, state);
}


AMREX_INLINE
void find_ignition_condition() {
    /*
      Main driver to determine the ignition condition
     */

    // Gravitational Acceleration for 1.4 solar mass and 11 km NS.
    // constexpr amrex::Real g = 1.5e14_rt;
    constexpr amrex::Real g = 2.45e14_rt;  // This is what settle version uses.
    constexpr amrex::Real R = 1.1e6_rt;

    State state;

    state.g = g;
    state.R = R;

    // Read in input parameters: X, Z, F_b, and mdot

    state.X0 = unit_test_rp::X;
    state.Z0 = unit_test_rp::Z;

    // Compute the Eddington Accretion Rate
    // mdot_Edd = 2 m_p c / [(1 + X) R σ_th ]
    // They used a value of 8.8e4 g cm^-2 s^-1
    // We can compute it at runtime in the future
    amrex::Real mdot_Edd = 8.8e4_rt;  // Come back later

    // Input mdot should be in eddington unit.
    state.mdot = unit_test_rp::mdot * mdot_Edd;

    // Input Fb is MeV/nucleon/mdot.
    // Brown & Bildstein showed that electron capture rates can
    // release ~ 1 MeV/nucleon deep in the crust
    // The compressional term can give additional c_p T ~ 20 T8 keV/nucleon
    // If don't include compressional term, use higher constant Fb
    // to account for compressional heating. They use ~1.5 MeV/nucleon
    // Convert to CGS.
    // Here they used mdot in Eddington unit
    state.Fb = unit_test_rp::Fb * unit_test_rp::mdot * 14.62_rt * 5.8e21_rt; // Not sure why multiplied by these constants, supposedly it is to convert to CGS...

    // Get the initial helium mass fraction
    state.Y0 = 1.0_rt - state.X0 - state.Z0;

    // Energy release per gram from hot CNO.
    // This is mass difference between helium and 4 hydrogen
    constexpr amrex::Real E_H = 6.4e18_rt;

    // Specific nuclear energy generation rate for hot CNO
    amrex::Real epsilon_H = 5.8e13 * (state.Z0 / 0.01);

    // Depletion Column Depth: yd = X mdot E_H / epsilon_H
    state.yd = state.X0 * state.mdot * E_H / epsilon_H;

    // Display input parameters to stdout
    std::cout << "Base flux is " << unit_test_rp::Fb << std::endl;
    std::cout << "In CGS, base flux is " << state.Fb << std::endl;
    std::cout << "Metallicity is " << unit_test_rp::Z << std::endl;
    std::cout << "Accreted H fraction is " << unit_test_rp::X << std::endl;
    std::cout << "mdot/Edd is " << unit_test_rp::mdot << std::endl;
    std::cout << std::scientific << "Depletion column depth is " << state.yd << std::endl;

    // Now compute the base column depth: yb via root find.
    // This is also the ignition column depth
    // Assume that we work with log10(yb)

    std::cout << "\n" << std::endl;
    std::cout << "Searching for ignition column depth..." << std::endl;
    find_yb(state);

    // Use solution to find density via EOS

    eos_t eos_state;
    eos_state.T = state.Tb;
    eos_state.p = state.g * state.yb;

    amrex::Real X = std::max(state.X0 * (1.0_rt - state.yb / state.yd), 0.0_rt);
    eos_state.xn[Species::H1-1] = X;
    eos_state.xn[Species::He4-1] = 1 - X - state.Z0;
    eos_state.xn[Species::O14-1] = O14 * state.Z0;
    eos_state.xn[Species::O15-1] = O15 * state.Z0;

    // Find density via Temperature and pressure
    eos(eos_input_tp, eos_state);

    std::cout << "\n" << std::endl;
    std::cout << "------------------------ SOLUTION ------------------------" << std::endl;
    std::cout << "Solved ignition column depth: y = " << state.yb << std::endl;
    std::cout << "With ignition temperature: T = " << state.Tb << std::endl;
    std::cout << "With ignition pressure: P = " << state.g * state.yb << std::endl;
    std::cout << "With ignition density: ϱ = " << eos_state.rho << std::endl;
    std::cout << "With relative height: z = " << state.z << std::endl;

    // Redshift factor used in settle
    constexpr amrex::Real ZZ = 1.31;

    std::cout << "With recurrence time (hours): y/mdot = " << ZZ * state.yb / state.mdot / 3600.0 << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "X(" << short_spec_names_cxx[n] << ") is " << eos_state.xn[n] << std::endl;
    }
}

#endif
