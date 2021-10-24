#ifndef rhs_utilities_H
#define rhs_utilities_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>

#include <network_utilities.H>
#include <actual_network.H>
#include <burn_type.H>
#include <screen.H>
#include <sneut5.H>
#include <jacobian_utilities.H>

// Forward declarations

namespace RHS
{

constexpr Real tab_tlo = 6.0e0_rt;
constexpr Real tab_thi = 10.0e0_rt;
constexpr int tab_per_decade = 2000;
constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr Real tab_tstp = (tab_thi - tab_tlo) / static_cast<Real>(tab_imax - 1);

extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRates, 1, 2, 1, nrattab> rattab;
extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRates, 1, 2, 1, nrattab> drattabdt;
extern AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;

struct rate_t
{
    // the forward and reverse rate data
    amrex::Array1D<amrex::Real, 1, Rates::NumRates> forward;
    amrex::Array1D<amrex::Real, 1, Rates::NumRates> reverse;
};

struct rate_tab_t
{
    Real alfa, beta, gama, delt;
    int iat;

    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void initialize (Real temp) {
        // hash locate
        constexpr int mp = 4;
        iat = static_cast<int>((std::log10(temp) - tab_tlo) / tab_tstp) + 1;
        iat = amrex::max(1, amrex::min(iat - 1, tab_imax - mp + 1));

        // setup the lagrange interpolation coefficients for a cubic
        Real x  = temp;
        Real x1 = ttab(iat);
        Real x2 = ttab(iat+1);
        Real x3 = ttab(iat+2);
        Real x4 = ttab(iat+3);
        Real a  = x - x1;
        Real b  = x - x2;
        Real c  = x - x3;
        Real d  = x - x4;
        Real e  = x1 - x2;
        Real f  = x1 - x3;
        Real g  = x1 - x4;
        Real h  = x2 - x3;
        Real p  = x2 - x4;
        Real q  = x3 - x4;
        alfa =  b * c * d / (e * f * g);
        beta = -a * c * d / (e * h * p);
        gama =  a * b * d / (f * h * q);
        delt = -a * b * c / (g * p * q);
    }
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE

// Internal implementation functions

namespace RHS_impl {

    // Calculate an integer factorial.
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    constexpr int factorial (int n)
    {
        if (n <= 1) {
            return 1;
        }
        else {
            return n * factorial(n - 1);
        }
    }

} // namespace RHS_impl

// Implementation of "constexpr for" based on
// https://artificial-mind.net/blog/2020/10/31/constexpr-for
//
// Approximates what one would get from a compile-time
// unrolling of the loop
// for (int i = 0; i < N; ++i) {
//    f(i);
// }
//
// The mechanism is recursive, we evaluate f(i) at the current
// i and then call the for loop at i+1. f() is a lambda function
// that provides the body of the loop and takes only an integer
// i as its argument. It is assumed that at the loop call site
// the lambda is declared [&] so that it can operate on whatever
// data it needs.

template<auto I, auto N, class F>
constexpr void constexpr_for (F&& f)
{
    if constexpr (I < N) {
        f(std::integral_constant<decltype(I), I>());
        constexpr_for<I+1, N>(f);
    }
}

// Determine if a rate is used in the RHS for a given species
// by seeing if its prefactor is nonzero.
template<int species, int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int is_rate_used ()
{
    constexpr rhs_t data = rhs_data(rate);

    static_assert(species >= 1 && species <= NumSpec);
    static_assert(rate >= 1 && rate <= Rates::NumRates);

    if (data.species_A == species ||
        data.species_B == species ||
        data.species_C == species ||
        data.species_D == species ||
        data.species_E == species ||
        data.species_F == species) {
        // Exclude intermediate rates that don't appear in the RHS.
        // We can identify these by the presence of an "extra" species
        // whose numerical ID is > NumSpec.
        if (data.species_A > NumSpec ||
            data.species_B > NumSpec ||
            data.species_C > NumSpec ||
            data.species_D > NumSpec ||
            data.species_E > NumSpec ||
            data.species_F > NumSpec) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return 0;
    }
}

// Calculate the density dependence term for tabulated rates. The RHS has a term
// that goes as rho**(exp_A + exp_B + exp_C) / rho (the denominator is because the
// LHS is X, not rho * X).
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int density_exponent_forward ()
{
    constexpr rhs_t data = rhs_data(rate);

    int exponent = 0;

    if constexpr (data.species_A >= 0) {
        exponent += data.exponent_A;
    }
    if constexpr (data.species_B >= 0) {
        exponent += data.exponent_B;
    }
    if constexpr (data.species_C >= 0) {
        exponent += data.exponent_C;
    }

    if (exponent > 0) {
        exponent -= 1;
    }

    return exponent;
}

// Same as the above but for the reverse reaction.
template <int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int density_exponent_reverse ()
{
    constexpr rhs_t data = rhs_data(rate);

    int exponent = 0;

    if constexpr (data.species_D >= 0) {
        exponent += data.exponent_D;
    }
    if constexpr (data.species_E >= 0) {
        exponent += data.exponent_E;
    }
    if constexpr (data.species_F >= 0) {
        exponent += data.exponent_F;
    }

    if (exponent > 0) {
        exponent -= 1;
    }

    return exponent;
}

template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_density_scaling (const Real& rho, Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    constexpr int forward_exponent = density_exponent_forward<rate>();
    constexpr int reverse_exponent = density_exponent_reverse<rate>();

    Real density_term_forward = std::pow(rho, forward_exponent);
    Real density_term_reverse = std::pow(rho, reverse_exponent);

    fr *= density_term_forward;
    frdt *= density_term_forward;
    rr *= density_term_reverse;
    rrdt *= density_term_reverse;
}

template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_screening (const plasma_state_t& pstate, Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    // The screening behavior depends on the type of reaction. We provide screening
    // here for the reaction classes we know about, and any other reactions are unscreened.

    constexpr rhs_t data = rhs_data(rate);

    if constexpr (data.number_A == 1 && data.number_B == 1 && data.number_C == 0) {
        // Forward reaction is A + B, screen using these two species

        constexpr Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr Real A1 = NetworkProperties::aion(data.species_A);

        constexpr Real Z2 = NetworkProperties::zion(data.species_B);
        constexpr Real A2 = NetworkProperties::aion(data.species_B);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        Real sc, scdt;
        actual_screen5(pstate, scn_fac, sc, scdt);

        frdt = frdt * sc + fr * scdt;
        fr   = fr * sc;
    }

    if constexpr (data.number_A == 2 && data.number_B == 0 && data.number_C == 0) {
        // Forward reaction is A + A, screen using just this species

        constexpr Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr Real A1 = NetworkProperties::aion(data.species_A);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        Real sc, scdt;
        actual_screen5(pstate, scn_fac, sc, scdt);

        frdt = frdt * sc + fr * scdt;
        fr   = fr * sc;
    }

    if constexpr (data.number_A == 3 && data.number_B == 0 && data.number_C == 0) {
        // Forward reaction is triple alpha or an equivalent, screen using A + A
        // and then A + X where X has twice the number of protons and neutrons.

        constexpr Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr Real A1 = NetworkProperties::aion(data.species_A);

        constexpr auto scn_fac1 = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        Real sc1, sc1dt;
        actual_screen5(pstate, scn_fac1, sc1, sc1dt);

        constexpr Real Z2 = 2.0_rt * Z1;
        constexpr Real A2 = 2.0_rt * A1;

        constexpr auto scn_fac2 = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        Real sc2, sc2dt;
        actual_screen5(pstate, scn_fac2, sc2, sc2dt);

        // Compute combined screening factor

        Real sc = sc1 * sc2;
        Real scdt = sc1dt * sc2 + sc1 * sc2dt;

        frdt = frdt * sc + fr * scdt;
        fr   = fr * sc;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void tabulate_rates ()
{
    using namespace Rates;

    for (int i = 1; i <= tab_imax; ++i) {
        Real temp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
        temp = std::pow(10.0e0_rt, temp);

        ttab(i) = temp;

        // Get the temperature factors
        tf_t tf = get_tfactors(temp);

        constexpr_for<1, NumRates+1>([&] (auto n)
        {
            constexpr int rate = n;

            Real fr = 0.0_rt;
            Real frdt = 0.0_rt;

            Real rr = 0.0_rt;
            Real rrdt = 0.0_rt;

            evaluate_analytical_rate<rate>(tf, fr, frdt, rr, rrdt);

            rattab(rate, 1, i)    = fr;
            rattab(rate, 2, i)    = rr;
            drattabdt(rate, 1, i) = frdt;
            drattabdt(rate, 2, i) = rrdt;
       });
    }
}

template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_tabulated_rate (const rate_tab_t& tab, Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    fr = (tab.alfa * rattab(rate, 1, tab.iat  ) +
          tab.beta * rattab(rate, 1, tab.iat+1) +
          tab.gama * rattab(rate, 1, tab.iat+2) +
          tab.delt * rattab(rate, 1, tab.iat+3));

    rr = (tab.alfa * rattab(rate, 2, tab.iat  ) +
          tab.beta * rattab(rate, 2, tab.iat+1) +
          tab.gama * rattab(rate, 2, tab.iat+2) +
          tab.delt * rattab(rate, 2, tab.iat+3));

    frdt = (tab.alfa * drattabdt(rate, 1, tab.iat  ) +
            tab.beta * drattabdt(rate, 1, tab.iat+1) +
            tab.gama * drattabdt(rate, 1, tab.iat+2) +
            tab.delt * drattabdt(rate, 1, tab.iat+3));

    rrdt = (tab.alfa * drattabdt(rate, 2, tab.iat  ) +
            tab.beta * drattabdt(rate, 2, tab.iat+1) +
            tab.gama * drattabdt(rate, 2, tab.iat+2) +
            tab.delt * drattabdt(rate, 2, tab.iat+3));
}

// Calculate the RHS term for a given species and rate.
//
// The general form of a reaction is
// n_a A + n_b B + n_c C <-> n_d D + n_e E + n_f F
// for species A, B, C, D, E, and F, where n_a particles of A,
// n_b particles of B, and n_C particles of C are consumed in
// the forward reaction and produced in the reverse reaction.
//
// For a given species, such as species A, the forward reaction
// term is of the form
// -n_A * Y(A)**a * Y(B)**b * Y(C)**c * forward_rate,
// and the reverse reaction term is of the form
//  n_A * Y(D)**d * Y(E)**e * Y(F)**f * reverse_rate.
// Here a, b, and c are reaction-specific exponents which usually,
// but not always, are equal to n_a, n_b, and n_c respectively.
//
// For example, in C12 + He4 <-> O16, species A is C12, species B
// is He4, species D is O16 (the other species are unused). Then
// n_a = n_b = n_d = 1, and a = b = d = 1. In the triple alpha forward
// reaction we have A = He4, D = C12, n_a = 3, a = 3, and n_d = 1.
//
// We assume the reaction rates do not include the identical particle
// factor, so we account for that here by dividing the rate by n!
// Note that we use the exponent to determine this factorial term, not
// the number consumed, because there are some reactions with special handling
// like the Si28 + 7 * He4 <-> Ni56 reaction in iso7, where the number
// of He4 consumed is not directly related to the actual reaction rate.
// In every other case than that, the exponent should be equal to the
// number consumed/produced, so there would be no difference.
//
// If a given reaction uses fewer than three species, we infer
// this by calling its index -1 and then not accessing it
// in the multiplication.
template<int species, int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real rhs_term (const burn_t& state, const Real& fr, const Real& rr)
{
    constexpr rhs_t data = rhs_data(rate);

    // First, compute the Y * rate component of both the forward and
    // reverse reactions, which is the same regardless of which species
    // we're producing or consuming.

    Real forward_term = fr;

    if constexpr (data.species_A >= 0) {
        Real Y_A = state.xn[data.species_A-1] * aion_inv[data.species_A-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(data.exponent_A);
        forward_term *= std::pow(Y_A, data.exponent_A) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (data.species_B >= 0) {
        Real Y_B = state.xn[data.species_B-1] * aion_inv[data.species_B-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(data.exponent_B);
        forward_term *= std::pow(Y_B, data.exponent_B) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (data.species_C >= 0) {
        Real Y_C = state.xn[data.species_C-1] * aion_inv[data.species_C-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(data.exponent_C);
        forward_term *= std::pow(Y_C, data.exponent_C) /
                        static_cast<Real>(identical_particle_factor);
    }

    Real reverse_term = rr;

    if constexpr (data.species_D >= 0) {
        Real Y_D = state.xn[data.species_D-1] * aion_inv[data.species_D-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(data.exponent_D);
        reverse_term *= std::pow(Y_D, data.exponent_D) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (data.species_E >= 0) {
        Real Y_E = state.xn[data.species_E-1] * aion_inv[data.species_E-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(data.exponent_E);
        reverse_term *= std::pow(Y_E, data.exponent_E) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (data.species_F >= 0) {
        Real Y_F = state.xn[data.species_F-1] * aion_inv[data.species_F-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(data.exponent_F);
        reverse_term *= std::pow(Y_F, data.exponent_F) /
                        static_cast<Real>(identical_particle_factor);
    }

    // Now compute the total contribution to this species.

    Real term = 0.0_rt;

    if constexpr (data.species_A == species) {
        term = data.number_A * (reverse_term - forward_term);
    }
    if constexpr (data.species_B == species) {
        term = data.number_B * (reverse_term - forward_term);
    }
    if constexpr (data.species_C == species) {
        term = data.number_C * (reverse_term - forward_term);
    }
    if constexpr (data.species_D == species) {
        term = data.number_D * (forward_term - reverse_term);
    }
    if constexpr (data.species_E == species) {
        term = data.number_E * (forward_term - reverse_term);
    }
    if constexpr (data.species_F == species) {
        term = data.number_F * (forward_term - reverse_term);
    }

    return term;
}

// Calculate the j'th Jacobian term for d(f(species1)) / d(species2).
//
// This follows the same scheme as the RHS. The forward term only
// has a contribution if the species we're taking the derivative
// with respect to (spec2) is one of (A, B, C). For the species
// that is spec2, we take the derivative by multiplying by the
// current exponent and then decrementing the exponent in the
// term. The same is done for the reverse term, which only has a
// contribution if spec2 is one of (D, E, F).
template<int spec1, int spec2, int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real jac_term (const burn_t& state, const Real& fr, const Real& rr)
{
    constexpr rhs_t data = rhs_data(rate);

    Real forward_term = 0.0_rt;

    if constexpr (spec2 == data.species_A || spec2 == data.species_B || spec2 == data.species_C) {

        forward_term = fr;

        if constexpr (data.species_A >= 0) {
            Real Y_A = state.xn[data.species_A-1] * aion_inv[data.species_A-1];

            constexpr int exponent = data.exponent_A;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            // Compute the forward term. It only has a contribution if
            // the species we're taking the derivative with respect to
            // (spec2) is one of (A, B, C). For the species that is spec2,
            // we take the derivative by multiplying by the current exponent
            // and then decrementing the exponent in the term.

            int exp = exponent;

            if constexpr (spec2 == data.species_A) {
                forward_term *= exp;
                --exp;
            }

            forward_term *= std::pow(Y_A, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (data.species_B >= 0) {
            Real Y_B = state.xn[data.species_B-1] * aion_inv[data.species_B-1];

            constexpr int exponent = data.exponent_B;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == data.species_B) {
                forward_term *= exp;
                --exp;
            }

            forward_term *= std::pow(Y_B, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (data.species_C >= 0) {
            Real Y_C = state.xn[data.species_C-1] * aion_inv[data.species_C-1];

            constexpr int exponent = data.exponent_C;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == data.species_C) {
                forward_term *= exp;
                --exp;
            }

            forward_term *= std::pow(Y_C, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

    }

    Real reverse_term = 0.0_rt;

    if constexpr (spec2 == data.species_D || spec2 == data.species_E || spec2 == data.species_F) {

        reverse_term = rr;

        if constexpr (data.species_D >= 0) {
            Real Y_D = state.xn[data.species_D-1] * aion_inv[data.species_D-1];

            constexpr int exponent = data.exponent_D;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == data.species_D) {
                reverse_term *= exp;
                --exp;
            }

            reverse_term *= std::pow(Y_D, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (data.species_E >= 0) {
            Real Y_E = state.xn[data.species_E-1] * aion_inv[data.species_E-1];

            constexpr int exponent = data.exponent_E;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == data.species_E) {
                reverse_term *= exp;
                --exp;
            }

            reverse_term *= std::pow(Y_E, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (data.species_F >= 0) {
            Real Y_F = state.xn[data.species_F-1] * aion_inv[data.species_F-1];

            constexpr int exponent = data.exponent_F;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == data.species_F) {
                reverse_term *= exp;
                --exp;
            }

            reverse_term *= std::pow(Y_F, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

    }

    // Now compute the total contribution to this species.

    Real term = 0.0_rt;

    if constexpr (data.species_A == spec1) {
        term = data.number_A * (reverse_term - forward_term);
    }
    if constexpr (data.species_B == spec1) {
        term = data.number_B * (reverse_term - forward_term);
    }
    if constexpr (data.species_C == spec1) {
        term = data.number_C * (reverse_term - forward_term);
    }
    if constexpr (data.species_D == spec1) {
        term = data.number_D * (forward_term - reverse_term);
    }
    if constexpr (data.species_E == spec1) {
        term = data.number_E * (forward_term - reverse_term);
    }
    if constexpr (data.species_F == spec1) {
        term = data.number_F * (forward_term - reverse_term);
    }

    return term;
}

template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void construct_rate (const Real& rho, const tf_t& tf, const rate_tab_t& tab,
                     const Array1D<Real, 1, NumSpec>& y, const plasma_state_t& pstate,
                     Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    using namespace Species;
    using namespace Rates;

    if (use_tables) {
        evaluate_tabulated_rate<rate>(tab, fr, frdt, rr, rrdt);
    }
    else {
        evaluate_analytical_rate<rate>(tf, fr, frdt, rr, rrdt);
    }

    // Set the density dependence

    apply_density_scaling<rate>(rho, fr, frdt, rr, rrdt);

    apply_screening<rate>(pstate, fr, frdt, rr, rrdt);
}

template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_rate (const Real& rho, const tf_t& tf, const rate_tab_t& tab,
               const Array1D<Real, 1, NumSpec>& y, const plasma_state_t& pstate,
               Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    using namespace Species;
    using namespace Rates;

    fr = 0.0;
    frdt = 0.0;
    rr = 0.0;
    rrdt = 0.0;

    constexpr rhs_t data = rhs_data(rate);

    // We can skip this step if this is an intermediate rate
    // (that is, involves at least one "extra" species).

    if constexpr (data.species_A > NumSpec ||
                  data.species_B > NumSpec ||
                  data.species_C > NumSpec ||
                  data.species_D > NumSpec ||
                  data.species_E > NumSpec ||
                  data.species_F > NumSpec) {
        return;
    }

    // Construct the (screened and density-scaled) rate.

    construct_rate<rate>(rho, tf, tab, y, pstate,
                         fr, frdt, rr, rrdt);

    // Construct all intermediate rates needed to augment this reaction.
    // To keep the problem bounded we assume that there are no more than
    // two intermediate reactions needed, however it is straightforward
    // to add more as needed (or to generalize this into a loop over all
    // rates if really needed). If any are present, postprocess the rate
    // using the additional terms.

    Real fr1, fr1dt, rr1, rr1dt;
    Real fr2, fr2dt, rr2, rr2dt;

    constexpr int rate1 = data.additional_reaction_1;
    constexpr int rate2 = data.additional_reaction_2;

    if constexpr (rate1 >= 0 || rate2 >= 0) {
        if constexpr (rate1 >= 0) {
            construct_rate<rate1>(rho, tf, tab, y, pstate,
                                  fr1, fr1dt, rr1, rr1dt);
        }

        if constexpr (rate2 >= 0) {
            construct_rate<rate1>(rho, tf, tab, y, pstate,
                                  fr2, fr2dt, rr2, rr2dt);
        }

        postprocess_rate<rate>(rho, tf, y,
                               fr, frdt, rr, rrdt,
                               fr1, fr1dt, rr1, rr1dt,
                               fr2, fr2dt, rr2, rr2dt);
    }

    // Handle any branching ratios.

    fr *= data.forward_branching_ratio;
    rr *= data.reverse_branching_ratio;

    frdt *= data.forward_branching_ratio;
    rrdt *= data.reverse_branching_ratio;
}

template<int spec>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ener_gener_rate (Real const& dydt)
{
    using namespace network;

    // This is basically e = m c**2
    return dydt * mion(spec) * C::Legacy::enuc_conv2;
}

AMREX_INLINE
void rhs_init ()
{
    rates_init();

    if (use_tables)
    {
        amrex::Print() << std::endl << " Initializing rate table" << std::endl;
        tabulate_rates();
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    // Convert X to Y.
    Array1D<Real, 1, NumSpec> y;
    for (int n = 1; n <= NumSpec; ++n) {
        y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    // Set up the state data, which is the same for all screening factors.
    plasma_state_t pstate;
    fill_plasma_state(pstate, state.T, state.rho, y);

    // Initialize the rate temperature term.
    tf_t tf = get_tfactors(state.T);
    rate_tab_t tab;
    if (use_tables) {
        tab.initialize(state.T);
    }

    // Initialize the RHS terms.
    for (int n = 1; n <= NumSpec; ++n) {
        ydot(n) = 0.0;
    }

    // Loop over all rates, and then loop over all species, and if the
    // rate affects that given species, add its contribution to the RHS.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n1)
    {
        constexpr int rate = n1;

        Real fr, frdt, rr, rrdt;
        get_rate<rate>(state.rho, tf, tab, y, pstate, fr, frdt, rr, rrdt);

        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int species = n2;

            if constexpr (is_rate_used<species, rate>()) {
                ydot(species) += rhs_term<species, rate>(state, fr, rr);
            }
        });
    });

    // Evaluate the neutrino cooling.
    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Compute the energy RHS term.
    ydot(net_ienuc) = -sneut;

    constexpr_for<1, NumSpec+1>([&] (auto n)
    {
        constexpr int species = n;

        ydot(net_ienuc) += ener_gener_rate<species>(ydot(species));
    });
}

// Analytical Jacobian
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void jac (burn_t& state, ArrayUtil::MathArray2D<1, neqs, 1, neqs>& jac)
{
    // Convert X to Y.
    Array1D<Real, 1, NumSpec> y;
    for (int n = 1; n <= NumSpec; ++n) {
        y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    // Set up the state data, which is the same for all screening factors.
    plasma_state_t pstate;
    fill_plasma_state(pstate, state.T, state.rho, y);

    // Initialize the rate temperature term.
    tf_t tf = get_tfactors(state.T);
    rate_tab_t tab;
    if (use_tables) {
        tab.initialize(state.T);
    }

    // Initialize the Jacobian terms.
    for (int i = 1; i <= neqs; ++i) {
        for (int j = 1; j <= neqs; ++j) {
            jac(i,j) = 0.0;
        }
    }

    // Loop over rates and compute Jacobian terms.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n1)
    {
        constexpr int rate = n1;

        Real fr, frdt, rr, rrdt;
        get_rate<rate>(state.rho, tf, tab, y, pstate, fr, frdt, rr, rrdt);

        // Species Jacobian elements with respect to other species.
        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int spec1 = n2;

            constexpr_for<1, NumSpec+1>([&] (auto n3)
            {
                constexpr int spec2 = n3;

                if constexpr (is_rate_used<spec1, rate>()) {
                    jac(spec1, spec2) += jac_term<spec1, spec2, rate>(state, fr, rr);
                }
            });
        });

        // Evaluate the Jacobian elements with respect to temperature.
        // We'll convert them from d/dT to d/de later.
        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int species = n2;

            if constexpr (is_rate_used<species, rate>()) {
                jac(species, net_ienuc) += rhs_term<species, rate>(state, frdt, rrdt);
            }
        });
    });

    // Evaluate the neutrino cooling.
    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar,
           sneut, dsneutdt, dsneutdd, snuda, snudz);

    // Energy generation rate Jacobian elements with respect to species

    constexpr_for<1, NumSpec+1>([&] (auto j)
    {
        constexpr int species = j;

        Real b1 = (-state.abar * state.abar * snuda + (NetworkProperties::zion(species) - state.zbar) * state.abar * snudz);
        jac(net_ienuc, species) = -b1;

        constexpr_for<1, NumSpec+1>([&] (auto i)
        {
            constexpr int s = i;

            jac(net_ienuc, species) += ener_gener_rate<s>(jac(s, species));
        });
    });

    // Convert from d/dT to d/de.
    jac(net_ienuc, net_ienuc) = -dsneutdt;

    constexpr_for<1, NumSpec+1>([&] (auto n)
    {
        constexpr int species = n;

        jac(n, net_ienuc) = temperature_to_energy_jacobian(state, jac(n, net_ienuc));
        jac(net_ienuc, net_ienuc) += ener_gener_rate<species>(jac(n, net_ienuc));
    });
}

} // namespace RHS

#endif
