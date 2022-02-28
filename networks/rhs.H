#ifndef rhs_H
#define rhs_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>

#include <ArrayUtilities.H>
#include <rhs_type.H>
#include <actual_network.H>
#include <burn_type.H>
#ifdef SCREENING
#include <screen.H>
#endif
#ifdef NEUTRINOS
#include <sneut5.H>
#endif
#include <jacobian_utilities.H>
#include <vode_type.H>

#ifdef NEW_NETWORK_IMPLEMENTATION

// Forward declarations

namespace RHS
{

// Rate tabulation data.
constexpr Real tab_tlo = 6.0e0_rt;
constexpr Real tab_thi = 10.0e0_rt;
constexpr int tab_per_decade = 2000;
constexpr int nrattab = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr int tab_imax = static_cast<int>(tab_thi - tab_tlo) * tab_per_decade + 1;
constexpr Real tab_tstp = (tab_thi - tab_tlo) / static_cast<Real>(tab_imax - 1);

extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRates, 1, 2, 1, nrattab> rattab;
extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRates, 1, 2, 1, nrattab> drattabdt;
extern AMREX_GPU_MANAGED Array1D<Real, 1, nrattab> ttab;

// A struct that contains the terms needed to evaluate a tabulated rate.
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
AMREX_GPU_HOST_DEVICE AMREX_INLINE
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

// Determine the index of a given intermediate reaction. We use the
// order of the original rate definitions
    
// Counts up the number of intermediate reactions. An intermediate
// reaction is defined as any reaction which contributes to the
// construction of some other reaction. Note that an intermediate
// reaction may directly contribute to the RHS itself in addition
// to being used in other reactions.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int num_intermediate_reactions ()
{
    int rate_is_intermediate[Rates::NumRates] = {0};

    // Loop through rates and increment the counter for any additional
    // reactions used by that rate. We do this as a first step to keep
    // the algorithm linear in the number of rates.
    for (int rate = 1; rate <= Rates::NumRates; ++rate) {
        rhs_t data = RHS::rhs_data(rate);

        if (data.additional_reaction_1 >= 1) {
            rate_is_intermediate[data.additional_reaction_1 - 1] += 1;
        }

        if (data.additional_reaction_2 >= 1) {
            rate_is_intermediate[data.additional_reaction_2 - 1] += 1;
        }

        if (data.additional_reaction_3 >= 1) {
            rate_is_intermediate[data.additional_reaction_3 - 1] += 1;
        }
    }

    // Now count up all intermediate rates.
    int count = 0;

    for (int rate = 1; rate <= Rates::NumRates; ++rate) {
        if (rate_is_intermediate[rate - 1] > 0) {
            ++count;
        }
    }

    return count;
}

// Locate an intermediate rate in the intermediate rates array.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int
locate_intermediate_rate_index (int intermediate_rate)
{
    constexpr int num_intermediate = num_intermediate_reactions();

    int indices[num_intermediate] = {0};

    // We loop through all rates and add their intermediate
    // rates to the temporary storage array, skipping ones
    // that have already been added.
    for (int rate = 1; rate <= Rates::NumRates; ++rate) {
        rhs_t data = RHS::rhs_data(rate);

        if (data.additional_reaction_1 >= 1) {
            for (int n = 1; n <= num_intermediate; ++n) {
                if (indices[n-1] == data.additional_reaction_1) {
                    // This rate has already been counted, so we're done.
                    break;
                }
                else if (indices[n-1] == 0) {
                    // This slot has not been filled yet, so opportunistically fill it.
                    indices[n-1] = data.additional_reaction_1;
                    break;
                }
            }
        }

        if (data.additional_reaction_2 >= 1) {
            for (int n = 1; n <= num_intermediate; ++n) {
                if (indices[n-1] == data.additional_reaction_2) {
                    // This rate has already been counted, so we're done.
                    break;
                }
                else if (indices[n-1] == 0) {
                    // This slot has not been filled yet, so opportunistically fill it.
                    indices[n-1] = data.additional_reaction_2;
                    break;
                }
            }
        }

        if (data.additional_reaction_3 >= 1) {
            for (int n = 1; n <= num_intermediate; ++n) {
                if (indices[n-1] == data.additional_reaction_3) {
                    // This rate has already been counted, so we're done.
                    break;
                }
                else if (indices[n-1] == 0) {
                    // This slot has not been filled yet, so opportunistically fill it.
                    indices[n-1] = data.additional_reaction_3;
                    break;
                }
            }
        }
    }

    // Finally, loop through this array and return the index corresponding
    // to the requested rate.

    for (int n = 1; n <= num_intermediate; ++n) {
        if (indices[n-1] == intermediate_rate) {
            return n;
        }
    }

    // If we did not match any intermediate rates, that means we are dealing with
    // a reaction that is only used once; return a negative number to indicate no match.

    return -1;
}

template<int n1, int n2>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int is_jacobian_term_used ()
{

#ifndef STRANG
    // currently SDC uses a different ordering of the elements from Strang,
    // so we need to generalize the logic below to allow it to be reused
    return 1;
#else
    // If either term is a non-species component, assume it is used.
    if (n1 > NumSpec || n2 > NumSpec) {
        return 1;
    }

    int term_is_used = 0;

    // Loop through all rates and see if any rate touches both species.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n)
    {
        constexpr int rate = n;
        constexpr int spec1 = n1;
        constexpr int spec2 = n2;

        constexpr rhs_t data = rhs_data(rate);

        int is_spec_1_used = 0;

        if (data.species_A == spec1 ||
            data.species_B == spec1 ||
            data.species_C == spec1 ||
            data.species_D == spec1 ||
            data.species_E == spec1 ||
            data.species_F == spec1) {
            // Exclude intermediate rates that don't appear in the RHS.
            // We can identify these by the presence of an "extra" species
            // whose numerical ID is > NumSpec.
            if (data.species_A <= NumSpec ||
                data.species_B <= NumSpec ||
                data.species_C <= NumSpec ||
                data.species_D <= NumSpec ||
                data.species_E <= NumSpec ||
                data.species_F <= NumSpec) {
                is_spec_1_used = 1;
            }
        }

        int is_spec_2_used = 0;

        if (data.species_A == spec2 ||
            data.species_B == spec2 ||
            data.species_C == spec2 ||
            data.species_D == spec2 ||
            data.species_E == spec2 ||
            data.species_F == spec2) {
            // Exclude intermediate rates that don't appear in the RHS.
            // We can identify these by the presence of an "extra" species
            // whose numerical ID is > NumSpec.
            if (data.species_A <= NumSpec ||
                data.species_B <= NumSpec ||
                data.species_C <= NumSpec ||
                data.species_D <= NumSpec ||
                data.species_E <= NumSpec ||
                data.species_F <= NumSpec) {
                is_spec_2_used = 1;
            }
        }

        if (is_spec_1_used || is_spec_2_used) {
            term_is_used = 1;
        }
    });

    return term_is_used;
#endif
}

typedef amrex::Array1D<int, 1, VODE_NEQS> IArray1D;
typedef amrex::Array1D<Real, 1, VODE_NEQS> RArray1D;
typedef ArrayUtil::MathArray2D<1, VODE_NEQS, 1, VODE_NEQS> RArray2D;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dgesl (RArray2D& a1, RArray1D& b1)
{
    auto const& a = reinterpret_cast<Array2D<Real, 0, VODE_NEQS-1, 0, VODE_NEQS-1>&>(a1);
    auto& b = reinterpret_cast<Array1D<Real, 0, VODE_NEQS-1>&>(b1);

    // solve a * x = b
    // first solve l * y = b
    constexpr_for<0, VODE_NEQS-1>([&] (auto n1)
    {
        constexpr int k = n1;

        Real t = b(k);
        constexpr_for<k+1, VODE_NEQS>([&] (auto n2)
        {
            constexpr int j = n2;

            if (is_jacobian_term_used<j, k>()) {
                b(j) += t * a(j,k);
            }
        });
    });

    // now solve u * x = y
    constexpr_for<0, VODE_NEQS>([&] (auto kb)
    {
        constexpr int k = VODE_NEQS - kb - 1;

        b(k) = b(k) / a(k,k);
        Real t = -b(k);

        constexpr_for<0, k>([&] (auto j)
        {
            if (is_jacobian_term_used<j, k>()) {
                b(j) += t * a(j,k);
            }
        });
    });
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dgefa (RArray2D& a1)
{
    auto& a = reinterpret_cast<Array2D<Real, 0, VODE_NEQS-1, 0, VODE_NEQS-1>&>(a1);

    // LU factorization in-place without pivoting.

    constexpr_for<0, VODE_NEQS-1>([&] (auto n1)
    {
        constexpr int k = n1;

        // compute multipliers
        Real t = -1.0_rt / a(k,k);
        constexpr_for<k+1, VODE_NEQS>([&] (auto n2)
        {
            constexpr int j = n2;

            if (is_jacobian_term_used<j, k>()) {
                a(j,k) *= t;
            }
        });

        // row elimination with column indexing
        constexpr_for<k+1, VODE_NEQS>([&] (auto n2)
        {
            constexpr int j = n2;

            t = a(k,j);
            constexpr_for<k+1, VODE_NEQS>([&] (auto n3)
            {
                constexpr int i = n3;

                if constexpr (is_jacobian_term_used<i, k>()) {
                    a(i,j) += t * a(i,k);
                }
            });
        });
    });
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

// Scale a rate using the density terms.
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_density_scaling (const Real& rho, Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    constexpr int forward_exponent = density_exponent_forward<rate>();
    constexpr int reverse_exponent = density_exponent_reverse<rate>();

    static_assert(forward_exponent <= 3);
    static_assert(reverse_exponent <= 3);

    // We know the exponent is an integer so we can construct the density term explicitly.
    Real density_term_forward = 1.0;
    Real density_term_reverse = 1.0;

    if constexpr (forward_exponent == 1) {
        density_term_forward = rho;
    }
    else if constexpr (forward_exponent == 2) {
        density_term_forward = rho * rho;
    }
    else if constexpr (forward_exponent == 3) {
        density_term_forward = rho * rho * rho;
    }

    if constexpr (reverse_exponent == 1) {
        density_term_reverse = rho;
    }
    else if constexpr (reverse_exponent == 2) {
        density_term_reverse = rho * rho;
    }
    else if constexpr (reverse_exponent == 3) {
        density_term_reverse = rho * rho * rho;
    }

    fr *= density_term_forward;
    frdt *= density_term_forward;
    rr *= density_term_reverse;
    rrdt *= density_term_reverse;
}

#ifdef SCREENING
// Apply the screening term to a given rate.
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_screening (const plasma_state_t& pstate, Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    // The screening behavior depends on the type of reaction. We provide screening
    // here for the reaction classes we know about, and any other reactions are unscreened.

    constexpr rhs_t data = rhs_data(rate);

    if constexpr (data.screen_forward_reaction == 0 && data.screen_reverse_reaction == 0) {
        return;
    }

    if constexpr (data.exponent_A == 1 && data.exponent_B == 1 && data.exponent_C == 0) {
        // Forward reaction is A + B, screen using these two species

        constexpr Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr Real A1 = NetworkProperties::aion(data.species_A);

        constexpr Real Z2 = NetworkProperties::zion(data.species_B);
        constexpr Real A2 = NetworkProperties::aion(data.species_B);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        // Insert a static assert (which will always pass) to require the
        // compiler to evaluate the screen factor at compile time.
        static_assert(scn_fac.z1 == Z1);

        Real sc, scdt;
        actual_screen5(pstate, scn_fac, sc, scdt);

        if constexpr (data.screen_forward_reaction == 1) {
            frdt = frdt * sc + fr * scdt;
            fr   = fr * sc;
        }

        if constexpr (data.screen_reverse_reaction == 1) {
            rrdt = rrdt * sc + rr * scdt;
            rr   = rr * sc;
        }
    }

    if constexpr (data.exponent_A == 2 && data.exponent_B == 0 && data.exponent_C == 0) {
        // Forward reaction is A + A, screen using just this species

        constexpr Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr Real A1 = NetworkProperties::aion(data.species_A);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        static_assert(scn_fac.z1 == Z1);

        Real sc, scdt;
        actual_screen5(pstate, scn_fac, sc, scdt);

        if constexpr (data.screen_forward_reaction == 1) {
            frdt = frdt * sc + fr * scdt;
            fr   = fr * sc;
        }

        if constexpr (data.screen_reverse_reaction == 1) {
            rrdt = rrdt * sc + rr * scdt;
            rr   = rr * sc;
        }
    }

    if constexpr (data.exponent_A == 3 && data.exponent_B == 0 && data.exponent_C == 0) {
        // Forward reaction is triple alpha or an equivalent, screen using A + A
        // and then A + X where X has twice the number of protons and neutrons.

        constexpr Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr Real A1 = NetworkProperties::aion(data.species_A);

        constexpr auto scn_fac1 = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        static_assert(scn_fac1.z1 == Z1);

        Real sc1, sc1dt;
        actual_screen5(pstate, scn_fac1, sc1, sc1dt);

        constexpr Real Z2 = 2.0_rt * Z1;
        constexpr Real A2 = 2.0_rt * A1;

        constexpr auto scn_fac2 = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac2.z1 == Z1);

        Real sc2, sc2dt;
        actual_screen5(pstate, scn_fac2, sc2, sc2dt);

        // Compute combined screening factor

        Real sc = sc1 * sc2;
        Real scdt = sc1dt * sc2 + sc1 * sc2dt;

        if constexpr (data.screen_forward_reaction == 1) {
            frdt = frdt * sc + fr * scdt;
            fr   = fr * sc;
        }

        if constexpr (data.screen_reverse_reaction == 1) {
            rrdt = rrdt * sc + rr * scdt;
            rr   = rr * sc;
        }
    }
}
#endif // SCREENING

// Apply the branching ratios to a given rate.
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_branching (Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    constexpr rhs_t data = rhs_data(rate);

    if constexpr (data.forward_branching_ratio != 1.0_rt) {
        fr *= data.forward_branching_ratio;
        frdt *= data.forward_branching_ratio;
    }

    if constexpr (data.reverse_branching_ratio != 1.0_rt) {
        rr *= data.reverse_branching_ratio;
        rrdt *= data.reverse_branching_ratio;
    }
}

// Do the initial tabulation of rates. We loop over various
// temperatures and evaluate all the rates at each temperature.
AMREX_INLINE
void tabulate_rates ()
{
    using namespace Rates;

    for (int i = 1; i <= tab_imax; ++i) {
        Real temp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
        temp = std::pow(10.0e0_rt, temp);

        ttab(i) = temp;

        // Get the temperature factors
        tf_t tf = get_tfactors(temp);

        // Arbitrary density, y, and y_e values (should be unused)
        const Real rho = 0.0_rt;
        const Real y_e = 0.0_rt;
        const Real eta = 0.0_rt;
        Array1D<Real, 1, NumSpec> y = {0.0_rt};

        constexpr_for<1, NumRates+1>([&] (auto n)
        {
            constexpr int rate = n;

            Real fr = 0.0_rt;
            Real frdt = 0.0_rt;

            Real rr = 0.0_rt;
            Real rrdt = 0.0_rt;

            constexpr rhs_t data = RHS::rhs_data(rate);

            if constexpr (data.rate_can_be_tabulated) {
                evaluate_analytical_rate<rate>(rho, tf, y_e, eta, y, fr, frdt, rr, rrdt);
            }

            rattab(rate, 1, i)    = fr;
            rattab(rate, 2, i)    = rr;
            drattabdt(rate, 1, i) = frdt;
            drattabdt(rate, 2, i) = rrdt;
       });
    }
}

// Evaluate a rate using the rate tables.
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_tabulated_rate (const Real& rho, const rate_tab_t& tab, const Real& y_e,
                              const Real& eta, const Array1D<Real, 1, NumSpec>& y,
                              Real& fr, Real& frdt, Real& rr, Real& rrdt)
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
// In most other cases than that, the exponent should be equal to the
// number consumed/produced, so there would be no difference. We allow
// a reaction to turn off the identical particle factor in case this scheme
// would cause problems; for example, in aprox19 the Fe52(a,p)Co55(g,p)Fe54
// reaction involves an equilibrium reaction so the identical particle factor
// does not apply despite the fact that two protons are involved.
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

        static_assert(data.exponent_A <= 3);

        Real df = 1.0;

        if constexpr (data.exponent_A == 1) {
            df = Y_A;
        }
        else if constexpr (data.exponent_A == 2) {
            df = Y_A * Y_A;
        }
        else if constexpr (data.exponent_A == 3) {
            df = Y_A * Y_A * Y_A;
        }

        if constexpr (data.apply_identical_particle_factor) {
            constexpr int identical_particle_factor = factorial(data.exponent_A);
            df *= 1.0_rt / static_cast<Real>(identical_particle_factor);
        }

        forward_term *= df;
    }

    if constexpr (data.species_B >= 0) {
        Real Y_B = state.xn[data.species_B-1] * aion_inv[data.species_B-1];

        static_assert(data.exponent_B <= 3);

        Real df = 1.0;

        if constexpr (data.exponent_B == 1) {
            df = Y_B;
        }
        else if constexpr (data.exponent_B == 2) {
            df = Y_B * Y_B;
        }
        else if constexpr (data.exponent_B == 3) {
            df = Y_B * Y_B * Y_B;
        }

        if constexpr (data.apply_identical_particle_factor) {
            constexpr int identical_particle_factor = factorial(data.exponent_B);
            df *= 1.0_rt / static_cast<Real>(identical_particle_factor);
        }

        forward_term *= df;
    }

    if constexpr (data.species_C >= 0) {
        Real Y_C = state.xn[data.species_C-1] * aion_inv[data.species_C-1];

        static_assert(data.exponent_C <= 3);

        Real df = 1.0;

        if constexpr (data.exponent_C == 1) {
            df = Y_C;
        }
        else if constexpr (data.exponent_C == 2) {
            df = Y_C * Y_C;
        }
        else if constexpr (data.exponent_C == 3) {
            df = Y_C * Y_C * Y_C;
        }

        if constexpr (data.apply_identical_particle_factor) {
            constexpr int identical_particle_factor = factorial(data.exponent_C);
            df *= 1.0_rt / static_cast<Real>(identical_particle_factor);
        }

        forward_term *= df;
    }

    Real reverse_term = rr;

    if constexpr (data.species_D >= 0) {
        Real Y_D = state.xn[data.species_D-1] * aion_inv[data.species_D-1];

        static_assert(data.exponent_D <= 3);

        Real dr = 1.0;

        if constexpr (data.exponent_D == 1) {
            dr = Y_D;
        }
        else if constexpr (data.exponent_D == 2) {
            dr = Y_D * Y_D;
        }
        else if constexpr (data.exponent_D == 3) {
            dr = Y_D * Y_D * Y_D;
        }

        if constexpr (data.apply_identical_particle_factor) {
            constexpr int identical_particle_factor = factorial(data.exponent_D);
            dr *= 1.0_rt / static_cast<Real>(identical_particle_factor);
        }

        reverse_term *= dr;
    }

    if constexpr (data.species_E >= 0) {
        Real Y_E = state.xn[data.species_E-1] * aion_inv[data.species_E-1];

        Real dr = 1.0;

        if constexpr (data.exponent_E == 1) {
            dr = Y_E;
        }
        else if constexpr (data.exponent_E == 2) {
            dr = Y_E * Y_E;
        }
        else if constexpr (data.exponent_E == 3) {
            dr = Y_E * Y_E * Y_E;
        }

        if constexpr (data.apply_identical_particle_factor) {
            constexpr int identical_particle_factor = factorial(data.exponent_E);
            dr *= 1.0_rt / static_cast<Real>(identical_particle_factor);
        }

        reverse_term *= dr;
    }

    if constexpr (data.species_F >= 0) {
        Real Y_F = state.xn[data.species_F-1] * aion_inv[data.species_F-1];

        Real dr = 1.0;

        if constexpr (data.exponent_F == 1) {
            dr = Y_F;
        }
        else if constexpr (data.exponent_F == 2) {
            dr = Y_F * Y_F;
        }
        else if constexpr (data.exponent_F == 3) {
            dr = Y_F * Y_F * Y_F;
        }

        if constexpr (data.apply_identical_particle_factor) {
            constexpr int identical_particle_factor = factorial(data.exponent_F);
            dr *= 1.0_rt / static_cast<Real>(identical_particle_factor);
        }

        reverse_term *= dr;
    }

    // Now compute the total contribution to this species.

    Real term = 0.0_rt;

    if constexpr (data.species_A == species) {
        term += data.number_A * (reverse_term - forward_term);
    }
    if constexpr (data.species_B == species) {
        term += data.number_B * (reverse_term - forward_term);
    }
    if constexpr (data.species_C == species) {
        term += data.number_C * (reverse_term - forward_term);
    }
    if constexpr (data.species_D == species) {
        term += data.number_D * (forward_term - reverse_term);
    }
    if constexpr (data.species_E == species) {
        term += data.number_E * (forward_term - reverse_term);
    }
    if constexpr (data.species_F == species) {
        term += data.number_F * (forward_term - reverse_term);
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

            // Compute the forward term. It only has a contribution if
            // the species we're taking the derivative with respect to
            // (spec2) is one of (A, B, C). For the species that is spec2,
            // we take the derivative by multiplying by the current exponent
            // and then decrementing the exponent in the term.

            constexpr int exp = (spec2 == data.species_A) ? exponent - 1 : exponent;

            static_assert(exp <= 3);

            Real df = 1.0;

            if constexpr (exp == 1) {
                df = Y_A;
            }
            else if constexpr (exp == 2) {
                df = Y_A * Y_A;
            }
            else if constexpr (exp == 3) {
                df = Y_A * Y_A * Y_A;
            }

            if constexpr (spec2 == data.species_A) {
                df *= exponent;
            }

            if constexpr (data.apply_identical_particle_factor) {
                constexpr int identical_particle_factor = factorial(exponent);
                df *= 1.0_rt / static_cast<Real>(identical_particle_factor);
            }

            forward_term *= df;
        }

        if constexpr (data.species_B >= 0) {
            Real Y_B = state.xn[data.species_B-1] * aion_inv[data.species_B-1];

            constexpr int exponent = data.exponent_B;

            constexpr int exp = (spec2 == data.species_B) ? exponent - 1 : exponent;

            static_assert(exp <= 3);

            Real df = 1.0;

            if constexpr (exp == 1) {
                df = Y_B;
            }
            else if constexpr (exp == 2) {
                df = Y_B * Y_B;
            }
            else if constexpr (exp == 3) {
                df = Y_B * Y_B * Y_B;
            }

            if constexpr (spec2 == data.species_B) {
                df *= exponent;
            }

            if constexpr (data.apply_identical_particle_factor) {
                constexpr int identical_particle_factor = factorial(exponent);
                df *= 1.0_rt / static_cast<Real>(identical_particle_factor);
            }

            forward_term *= df;
        }

        if constexpr (data.species_C >= 0) {
            Real Y_C = state.xn[data.species_C-1] * aion_inv[data.species_C-1];

            constexpr int exponent = data.exponent_C;

            constexpr int exp = (spec2 == data.species_C) ? exponent - 1 : exponent;

            static_assert(exp <= 3);

            Real df = 1.0;

            if constexpr (exp == 1) {
                df = Y_C;
            }
            else if constexpr (exp == 2) {
                df = Y_C * Y_C;
            }
            else if constexpr (exp == 3) {
                df = Y_C * Y_C * Y_C;
            }

            if constexpr (spec2 == data.species_C) {
                df *= exponent;
            }

            if constexpr (data.apply_identical_particle_factor) {
                constexpr int identical_particle_factor = factorial(exponent);
                df *= 1.0_rt / static_cast<Real>(identical_particle_factor);
            }

            forward_term *= df;
        }

    }

    Real reverse_term = 0.0_rt;

    if constexpr (spec2 == data.species_D || spec2 == data.species_E || spec2 == data.species_F) {

        reverse_term = rr;

        if constexpr (data.species_D >= 0) {
            Real Y_D = state.xn[data.species_D-1] * aion_inv[data.species_D-1];

            constexpr int exponent = data.exponent_D;

            constexpr int exp = (spec2 == data.species_D) ? exponent - 1 : exponent;

            static_assert(exp <= 3);

            Real dr = 1.0;

            if constexpr (exp == 1) {
                dr = Y_D;
            }
            else if constexpr (exp == 2) {
                dr = Y_D * Y_D;
            }
            else if constexpr (exp == 3) {
                dr = Y_D * Y_D * Y_D;
            }

            if constexpr (spec2 == data.species_D) {
                dr *= exponent;
            }

            if constexpr (data.apply_identical_particle_factor) {
                constexpr int identical_particle_factor = factorial(exponent);
                dr *= 1.0_rt / static_cast<Real>(identical_particle_factor);
            }

            reverse_term *= dr;
        }

        if constexpr (data.species_E >= 0) {
            Real Y_E = state.xn[data.species_E-1] * aion_inv[data.species_E-1];

            constexpr int exponent = data.exponent_E;

            constexpr int exp = (spec2 == data.species_E) ? exponent - 1 : exponent;

            static_assert(exp <= 3);

            Real dr = 1.0;

            if constexpr (exp == 1) {
                dr = Y_E;
            }
            else if constexpr (exp == 2) {
                dr = Y_E * Y_E;
            }
            else if constexpr (exp == 3) {
                dr = Y_E * Y_E * Y_E;
            }

            if constexpr (spec2 == data.species_E) {
                dr *= exponent;
            }

            if constexpr (data.apply_identical_particle_factor) {
                constexpr int identical_particle_factor = factorial(exponent);
                dr *= 1.0_rt / static_cast<Real>(identical_particle_factor);
            }

            reverse_term *= dr;
        }

        if constexpr (data.species_F >= 0) {
            Real Y_F = state.xn[data.species_F-1] * aion_inv[data.species_F-1];

            constexpr int exponent = data.exponent_F;

            constexpr int exp = (spec2 == data.species_F) ? exponent - 1 : exponent;

            static_assert(exp <= 3);

            Real dr = 1.0;

            if constexpr (exp == 1) {
                dr = Y_F;
            }
            else if constexpr (exp == 2) {
                dr = Y_F * Y_F;
            }
            else if constexpr (exp == 3) {
                dr = Y_F * Y_F * Y_F;
            }

            if constexpr (spec2 == data.species_F) {
                dr *= exponent;
            }

            if constexpr (data.apply_identical_particle_factor) {
                constexpr int identical_particle_factor = factorial(exponent);
                dr *= 1.0_rt / static_cast<Real>(identical_particle_factor);
            }

            reverse_term *= dr;
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
void construct_rate (const Real& rho, const tf_t& tf, const rate_tab_t& tab, const Real& y_e,
                     const Real& eta, const Array1D<Real, 1, NumSpec>& y, const plasma_state_t& pstate,
                     Real& fr, Real& frdt, Real& rr, Real& rrdt)
{
    using namespace Species;
    using namespace Rates;

    fr = 0.0;
    frdt = 0.0;
    rr = 0.0;
    rrdt = 0.0;

    constexpr rhs_t data = RHS::rhs_data(rate);

    if (use_tables && data.rate_can_be_tabulated) {
        evaluate_tabulated_rate<rate>(rho, tab, y_e, eta, y, fr, frdt, rr, rrdt);
    }
    else {
        evaluate_analytical_rate<rate>(rho, tf, y_e, eta, y, fr, frdt, rr, rrdt);
    }

    // Set the density dependence

    apply_density_scaling<rate>(rho, fr, frdt, rr, rrdt);

#ifdef SCREENING
    // Screen

    apply_screening<rate>(pstate, fr, frdt, rr, rrdt);
#endif

    // Branching ratios

    apply_branching<rate>(fr, frdt, rr, rrdt);
}

template<int rate, class Arr>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fill_additional_rates (Arr intermediate_fr, Arr intermediate_frdt,
                            Arr intermediate_rr, Arr intermediate_rrdt,
                            Real& fr1, Real& fr1dt, Real& rr1, Real& rr1dt,
                            Real& fr2, Real& fr2dt, Real& rr2, Real& rr2dt,
                            Real& fr3, Real& fr3dt, Real& rr3, Real& rr3dt)
{
    constexpr rhs_t data = RHS::rhs_data(rate);

    constexpr int rate1 = data.additional_reaction_1;
    constexpr int rate2 = data.additional_reaction_2;
    constexpr int rate3 = data.additional_reaction_3;

    if constexpr (rate1 >= 0) {
        constexpr int index = locate_intermediate_rate_index(rate1);

        static_assert(index >= 1 && index <= num_intermediate_reactions());

        fr1 = intermediate_fr(index);
        fr1dt = intermediate_frdt(index);
        rr1 = intermediate_rr(index);
        rr1dt = intermediate_rrdt(index);
    }

    if constexpr (rate2 >= 0) {
        constexpr int index = locate_intermediate_rate_index(rate2);

        static_assert(index >= 1 && index <= num_intermediate_reactions());

        fr2 = intermediate_fr(index);
        fr2dt = intermediate_frdt(index);
        rr2 = intermediate_rr(index);
        rr2dt = intermediate_rrdt(index);
    }

    if constexpr (rate3 >= 0) {
        constexpr int index = locate_intermediate_rate_index(rate3);

        static_assert(index >= 1 && index <= num_intermediate_reactions());

        fr3 = intermediate_fr(index);
        fr3dt = intermediate_frdt(index);
        rr3 = intermediate_rr(index);
        rr3dt = intermediate_rrdt(index);
    }
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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
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

    // Count up number of intermediate rates (rates that are used in any other reaction).
    constexpr int num_intermediate = num_intermediate_reactions();

    // We cannot have a zero-sized array, so just set the array size to 1 in that case.
    constexpr int intermediate_array_size = num_intermediate > 0 ? num_intermediate : 1;

    // Define forward and reverse (and d/dT) rate arrays.
    Array1D<Real, 1, intermediate_array_size> intermediate_fr;
    Array1D<Real, 1, intermediate_array_size> intermediate_frdt;
    Array1D<Real, 1, intermediate_array_size> intermediate_rr;
    Array1D<Real, 1, intermediate_array_size> intermediate_rrdt;

    // Fill all intermediate rates first.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n)
    {
        constexpr int rate = n;

        constexpr int index = locate_intermediate_rate_index(rate);

        if constexpr (index >= 1) {
            construct_rate<rate>(state.rho, tf, tab, state.y_e, state.eta, y, pstate,
                                 intermediate_fr(index),
                                 intermediate_frdt(index),
                                 intermediate_rr(index),
                                 intermediate_rrdt(index));
         }
    });

    // Loop over all rates, and then loop over all species, and if the
    // rate affects that given species, add its contribution to the RHS.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n1)
    {
        constexpr int rate = n1;

        Real fr, frdt, rr, rrdt;

        // We only need to compute the rate at this point if it's not intermediate. If it
        // is intermediate, retrieve it from the cached array.

        constexpr int index = locate_intermediate_rate_index(rate);
        if constexpr (index < 0) {
            construct_rate<rate>(state.rho, tf, tab, state.y_e, state.eta, y, pstate, fr, frdt, rr, rrdt);
        }
        else {
            fr = intermediate_fr(index);
            frdt = intermediate_frdt(index);
            rr = intermediate_rr(index);
            rrdt = intermediate_rrdt(index);
        }

        // Locate all intermediate rates needed to augment this reaction.
        // To keep the problem bounded we assume that there are no more than
        // three intermediate reactions needed.

        Real fr1, fr1dt, rr1, rr1dt;
        Real fr2, fr2dt, rr2, rr2dt;
        Real fr3, fr3dt, rr3, rr3dt;

        fill_additional_rates<rate>(intermediate_fr, intermediate_frdt,
                                    intermediate_rr, intermediate_rrdt,
                                    fr1, fr1dt, rr1, rr1dt,
                                    fr2, fr2dt, rr2, rr2dt,
                                    fr3, fr3dt, rr3, rr3dt);

        // Perform rate postprocessing, using additional reactions as inputs.
        // If there is no postprocessing for this rate, this will be a no-op.

        postprocess_rate<rate>(state.rho, tf, state.y_e, state.eta, y,
                               fr, frdt, rr, rrdt,
                               fr1, fr1dt, rr1, rr1dt,
                               fr2, fr2dt, rr2, rr2dt,
                               fr3, fr3dt, rr3, rr3dt);

        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int species = n2;

            if constexpr (is_rate_used<species, rate>()) {
                ydot(species) += rhs_term<species, rate>(state, fr, rr);
            }
        });
    });

    // Evaluate the neutrino cooling.
#ifdef NEUTRINOS
    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);
#else
    Real sneut = 0.0;
#endif

    // Compute the energy RHS term.
    ydot(net_ienuc) = -sneut;

    constexpr_for<1, NumSpec+1>([&] (auto n)
    {
        constexpr int species = n;

        ydot(net_ienuc) += ener_gener_rate<species>(ydot(species));
    });
}

// Analytical Jacobian
AMREX_GPU_HOST_DEVICE AMREX_INLINE
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

    // Count up number of intermediate rates (rates that are used in any other reaction).
    constexpr int num_intermediate = num_intermediate_reactions();

    // We cannot have a zero-sized array, so just set the array size to 1 in that case.
    constexpr int intermediate_array_size = num_intermediate > 0 ? num_intermediate : 1;

    // Define forward and reverse (and d/dT) rate arrays.
    Array1D<Real, 1, intermediate_array_size> intermediate_fr;
    Array1D<Real, 1, intermediate_array_size> intermediate_frdt;
    Array1D<Real, 1, intermediate_array_size> intermediate_rr;
    Array1D<Real, 1, intermediate_array_size> intermediate_rrdt;

    Real fr1, fr1dt, rr1, rr1dt;
    Real fr2, fr2dt, rr2, rr2dt;
    Real fr3, fr3dt, rr3, rr3dt;

    // Fill all intermediate rates first.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n)
    {
        constexpr int rate = n;

        constexpr int index = locate_intermediate_rate_index(rate);

        if constexpr (index >= 1) {
            construct_rate<rate>(state.rho, tf, tab, state.y_e, state.eta, y, pstate,
                                 intermediate_fr(index),
                                 intermediate_frdt(index),
                                 intermediate_rr(index),
                                 intermediate_rrdt(index));
         }
    });

    // Loop over rates and compute Jacobian terms.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n1)
    {
        constexpr int rate = n1;

        Real fr, frdt, rr, rrdt;

        // We only need to compute the rate at this point if it's not intermediate. If it
        // is intermediate, retrieve it from the cached array.

        constexpr int index = locate_intermediate_rate_index(rate);
        if constexpr (index < 0) {
            construct_rate<rate>(state.rho, tf, tab, state.y_e, state.eta, y, pstate, fr, frdt, rr, rrdt);
        }
        else {
            fr = intermediate_fr(index);
            frdt = intermediate_frdt(index);
            rr = intermediate_rr(index);
            rrdt = intermediate_rrdt(index);
        }

        // Locate all intermediate rates needed to augment this reaction.
        // To keep the problem bounded we assume that there are no more than
        // three intermediate reactions needed.

        fill_additional_rates<rate>(intermediate_fr, intermediate_frdt,
                                    intermediate_rr, intermediate_rrdt,
                                    fr1, fr1dt, rr1, rr1dt,
                                    fr2, fr2dt, rr2, rr2dt,
                                    fr3, fr3dt, rr3, rr3dt);

        // Perform rate postprocessing, using additional reactions as inputs.
        // If there is no postprocessing for this rate, this will be a no-op.

        postprocess_rate<rate>(state.rho, tf, state.y_e, state.eta, y,
                               fr, frdt, rr, rrdt,
                               fr1, fr1dt, rr1, rr1dt,
                               fr2, fr2dt, rr2, rr2dt,
                               fr3, fr3dt, rr3, rr3dt);

        // Species Jacobian elements with respect to other species.
        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int spec1 = n2;

            constexpr_for<1, NumSpec+1>([&] (auto n3)
            {
                if constexpr (is_rate_used<spec1, rate>()) {
                    constexpr int spec2 = n3;
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
#ifdef NEUTRINOS
    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar,
           sneut, dsneutdt, dsneutdd, snuda, snudz);
#else
    Real sneut = 0.0, dsneutdt = 0.0, dsneutdd = 0.0, snuda = 0.0, snudz = 0.0;
#endif

    jac(net_ienuc, net_ienuc) = -temperature_to_energy_jacobian(state, dsneutdt);

    constexpr_for<1, NumSpec+1>([&] (auto j)
    {
        constexpr int species = j;

        // Energy generation rate Jacobian elements with respect to species.
        Real b1 = (-state.abar * state.abar * snuda + (NetworkProperties::zion(species) - state.zbar) * state.abar * snudz);
        jac(net_ienuc, species) = -b1;

        constexpr_for<1, NumSpec+1>([&] (auto i)
        {
            constexpr int s = i;

            jac(net_ienuc, species) += ener_gener_rate<s>(jac(s, species));
        });

        // Convert previously computed terms from d/dT to d/de.
        jac(species, net_ienuc) = temperature_to_energy_jacobian(state, jac(species, net_ienuc));

        // Compute df(e) / de term.
        jac(net_ienuc, net_ienuc) += ener_gener_rate<species>(jac(species, net_ienuc));
    });
}

} // namespace RHS

// For legacy reasons, implement actual_rhs() and actual_jac() interfaces outside the RHS
// namespace. This should be retired later when nothing still depends on those names.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    RHS::rhs(state, ydot);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac (burn_t& state, ArrayUtil::MathArray2D<1, neqs, 1, neqs>& jac)
{
    RHS::jac(state, jac);
}

#endif // NEW_NETWORK_IMPLEMENTATION

#endif
