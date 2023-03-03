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
#include <integrator_data.H>

#ifdef NEW_NETWORK_IMPLEMENTATION

// Forward declarations

namespace RHS
{

// Rate tabulation data.
extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRates, 1, 2, 1, nrattab> rattab;
extern AMREX_GPU_MANAGED Array3D<Real, 1, Rates::NumRates, 1, 2, 1, nrattab> drattabdt;

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

    // num_intermediate may be 0, so use amrex::Array instead of a C-style array
    amrex::Array<int, num_intermediate> indices = {};

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

typedef amrex::Array1D<int, 1, INT_NEQS> IArray1D;
typedef amrex::Array1D<Real, 1, INT_NEQS> RArray1D;
typedef ArrayUtil::MathArray2D<1, INT_NEQS, 1, INT_NEQS> RArray2D;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void dgesl (RArray2D& a1, RArray1D& b1)
{
    auto const& a = reinterpret_cast<Array2D<Real, 0, INT_NEQS-1, 0, INT_NEQS-1>&>(a1);
    auto& b = reinterpret_cast<Array1D<Real, 0, INT_NEQS-1>&>(b1);

    // solve a * x = b
    // first solve l * y = b
    constexpr_for<0, INT_NEQS-1>([&] (auto n1)
    {
        constexpr int k = n1;

        Real t = b(k);
        constexpr_for<k+1, INT_NEQS>([&] (auto n2)
        {
            constexpr int j = n2;

            if (is_jacobian_term_used<j, k>()) {
                b(j) += t * a(j,k);
            }
        });
    });

    // now solve u * x = y
    constexpr_for<0, INT_NEQS>([&] (auto kb)
    {
        constexpr int k = INT_NEQS - kb - 1;

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
    auto& a = reinterpret_cast<Array2D<Real, 0, INT_NEQS-1, 0, INT_NEQS-1>&>(a1);

    // LU factorization in-place without pivoting.

    constexpr_for<0, INT_NEQS-1>([&] (auto n1)
    {
        [[maybe_unused]] constexpr int k = n1;

        // compute multipliers
        Real t = -1.0_rt / a(k,k);
        constexpr_for<k+1, INT_NEQS>([&] (auto n2)
        {
            [[maybe_unused]] constexpr int j = n2;

            if (is_jacobian_term_used<j, k>()) {
                a(j,k) *= t;
            }
        });

        // row elimination with column indexing
        constexpr_for<k+1, INT_NEQS>([&] (auto n2)
        {
            [[maybe_unused]] constexpr int j = n2;

            t = a(k,j);
            constexpr_for<k+1, INT_NEQS>([&] (auto n3)
            {
                [[maybe_unused]] constexpr int i = n3;

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
void apply_density_scaling (const rhs_state_t& state, rate_t& rates)
{
    constexpr int forward_exponent = density_exponent_forward<rate>();
    constexpr int reverse_exponent = density_exponent_reverse<rate>();

    static_assert(forward_exponent <= 3);
    static_assert(reverse_exponent <= 3);

    // We know the exponent is an integer so we can construct the density term explicitly.
    Real density_term_forward = 1.0;
    Real density_term_reverse = 1.0;

    if constexpr (forward_exponent == 1) {
        density_term_forward = state.rho;
    }
    else if constexpr (forward_exponent == 2) {
        density_term_forward = state.rho * state.rho;
    }
    else if constexpr (forward_exponent == 3) {
        density_term_forward = state.rho * state.rho * state.rho;
    }

    if constexpr (reverse_exponent == 1) {
        density_term_reverse = state.rho;
    }
    else if constexpr (reverse_exponent == 2) {
        density_term_reverse = state.rho * state.rho;
    }
    else if constexpr (reverse_exponent == 3) {
        density_term_reverse = state.rho * state.rho * state.rho;
    }

    rates.fr *= density_term_forward;
    rates.frdt *= density_term_forward;
    rates.rr *= density_term_reverse;
    rates.rrdt *= density_term_reverse;
}

#ifdef SCREENING
// Apply the screening term to a given rate.
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_screening (const rhs_state_t& state, rate_t& rates)
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
        actual_screen(state.pstate, scn_fac, sc, scdt);

        if constexpr (data.screen_forward_reaction == 1) {
            rates.frdt = rates.frdt * sc + rates.fr * scdt;
            rates.fr   = rates.fr * sc;
        }

        if constexpr (data.screen_reverse_reaction == 1) {
            rates.rrdt = rates.rrdt * sc + rates.rr * scdt;
            rates.rr   = rates.rr * sc;
        }
    }

    if constexpr (data.exponent_A == 2 && data.exponent_B == 0 && data.exponent_C == 0) {
        // Forward reaction is A + A, screen using just this species

        constexpr Real Z1 = NetworkProperties::zion(data.species_A);
        constexpr Real A1 = NetworkProperties::aion(data.species_A);

        constexpr auto scn_fac = scrn::calculate_screen_factor(Z1, A1, Z1, A1);

        static_assert(scn_fac.z1 == Z1);

        Real sc, scdt;
        actual_screen(state.pstate, scn_fac, sc, scdt);

        if constexpr (data.screen_forward_reaction == 1) {
            rates.frdt = rates.frdt * sc + rates.fr * scdt;
            rates.fr   = rates.fr * sc;
        }

        if constexpr (data.screen_reverse_reaction == 1) {
            rates.rrdt = rates.rrdt * sc + rates.rr * scdt;
            rates.rr   = rates.rr * sc;
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
        actual_screen(state.pstate, scn_fac1, sc1, sc1dt);

        constexpr Real Z2 = 2.0_rt * Z1;
        constexpr Real A2 = 2.0_rt * A1;

        constexpr auto scn_fac2 = scrn::calculate_screen_factor(Z1, A1, Z2, A2);

        static_assert(scn_fac2.z1 == Z1);

        Real sc2, sc2dt;
        actual_screen(state.pstate, scn_fac2, sc2, sc2dt);

        // Compute combined screening factor

        Real sc = sc1 * sc2;
        Real scdt = sc1dt * sc2 + sc1 * sc2dt;

        if constexpr (data.screen_forward_reaction == 1) {
            rates.frdt = rates.frdt * sc + rates.fr * scdt;
            rates.fr   = rates.fr * sc;
        }

        if constexpr (data.screen_reverse_reaction == 1) {
            rates.rrdt = rates.rrdt * sc + rates.rr * scdt;
            rates.rr   = rates.rr * sc;
        }
    }
}
#endif // SCREENING

// Apply the branching ratios to a given rate.
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_branching (rate_t& rates)
{
    constexpr rhs_t data = rhs_data(rate);

    if constexpr (data.forward_branching_ratio != 1.0_rt) {
        rates.fr *= data.forward_branching_ratio;
        rates.frdt *= data.forward_branching_ratio;
    }

    if constexpr (data.reverse_branching_ratio != 1.0_rt) {
        rates.rr *= data.reverse_branching_ratio;
        rates.rrdt *= data.reverse_branching_ratio;
    }
}

// Do the initial tabulation of rates. We loop over various
// temperatures and evaluate all the rates at each temperature.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void tabulate_rates ()
{
    using namespace Rates;

    for (int i = 1; i <= tab_imax; ++i) {
        Real temp = tab_tlo + static_cast<Real>(i-1) * tab_tstp;
        temp = std::pow(10.0e0_rt, temp);

        ttab(i) = temp;

        rhs_state_t state;

        // Get the temperature factors
        state.tf = get_tfactors(temp);

        // Arbitrary density, y, and y_e values (should be unused)
        state.rho = 0.0_rt;
        state.y_e = 0.0_rt;
        state.eta = 0.0_rt;
        state. y = {0.0_rt};

        constexpr_for<1, NumRates+1>([&] (auto n)
        {
            [[maybe_unused]] constexpr int rate = n;

            rate_t rates;
            rates.fr = 0.0_rt;
            rates.rr = 0.0_rt;
            rates.frdt = 0.0_rt;
            rates.rrdt = 0.0_rt;

            constexpr rhs_t data = RHS::rhs_data(rate);

            if constexpr (data.rate_can_be_tabulated) {
                evaluate_analytical_rate<rate>(state, rates);
            }

            rattab(rate, 1, i)    = rates.fr;
            rattab(rate, 2, i)    = rates.rr;
            drattabdt(rate, 1, i) = rates.frdt;
            drattabdt(rate, 2, i) = rates.rrdt;
       });
    }
}

// Evaluate a rate using the rate tables.
template<int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_tabulated_rate (const rhs_state_t& state, rate_t& rates)
{
    rates.fr = (state.tab.alfa * rattab(rate, 1, state.tab.iat  ) +
                state.tab.beta * rattab(rate, 1, state.tab.iat+1) +
                state.tab.gama * rattab(rate, 1, state.tab.iat+2) +
                state.tab.delt * rattab(rate, 1, state.tab.iat+3));

    rates.rr = (state.tab.alfa * rattab(rate, 2, state.tab.iat  ) +
                state.tab.beta * rattab(rate, 2, state.tab.iat+1) +
                state.tab.gama * rattab(rate, 2, state.tab.iat+2) +
                state.tab.delt * rattab(rate, 2, state.tab.iat+3));

    rates.frdt = (state.tab.alfa * drattabdt(rate, 1, state.tab.iat  ) +
                  state.tab.beta * drattabdt(rate, 1, state.tab.iat+1) +
                  state.tab.gama * drattabdt(rate, 1, state.tab.iat+2) +
                  state.tab.delt * drattabdt(rate, 1, state.tab.iat+3));

    rates.rrdt = (state.tab.alfa * drattabdt(rate, 2, state.tab.iat  ) +
                  state.tab.beta * drattabdt(rate, 2, state.tab.iat+1) +
                  state.tab.gama * drattabdt(rate, 2, state.tab.iat+2) +
                  state.tab.delt * drattabdt(rate, 2, state.tab.iat+3));
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
template<int species, int rate, int use_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr std::pair<Real, Real> rhs_term (const burn_t& state, const rate_t& rates)
{
    constexpr rhs_t data = rhs_data(rate);

    // First, compute the Y * rate component of both the forward and
    // reverse reactions, which is the same regardless of which species
    // we're producing or consuming.

    Real forward_term = use_T_derivatives ? rates.frdt : rates.fr;

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

    Real reverse_term = use_T_derivatives ? rates.rrdt : rates.rr;

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

    if constexpr (data.species_A == species) {
        forward_term *= -data.number_A;
        reverse_term *= data.number_A;
    }
    else if constexpr (data.species_B == species) {
        forward_term *= -data.number_B;
        reverse_term *= data.number_B;
    }
    else if constexpr (data.species_C == species) {
        forward_term *= -data.number_C;
        reverse_term *= data.number_C;
    }
    else if constexpr (data.species_D == species) {
        forward_term *= data.number_D;
        reverse_term *= -data.number_D;
    }
    else if constexpr (data.species_E == species) {
        forward_term *= data.number_E;
        reverse_term *= -data.number_E;
    }
    else if constexpr (data.species_F == species) {
        forward_term *= data.number_F;
        reverse_term *= -data.number_F;
    }
    else {
        forward_term = 0.0;
        reverse_term = 0.0;
    }

    return {forward_term, reverse_term};
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
constexpr Real jac_term (const burn_t& state, const rate_t& rates)
{
    constexpr rhs_t data = rhs_data(rate);

    Real forward_term = 0.0_rt;

    if constexpr (is_rate_used<spec1, rate>() &&
                  (spec2 == data.species_A || spec2 == data.species_B || spec2 == data.species_C)) {

        forward_term = rates.fr;

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

    if constexpr (is_rate_used<spec1, rate>() &&
                  (spec2 == data.species_D || spec2 == data.species_E || spec2 == data.species_F)) {

        reverse_term = rates.rr;

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
void construct_rate (const rhs_state_t& state, rate_t& rates)
{
    using namespace Species;
    using namespace Rates;

    rates.fr = 0.0;
    rates.frdt = 0.0;
    rates.rr = 0.0;
    rates.rrdt = 0.0;

    constexpr rhs_t data = RHS::rhs_data(rate);

    if (use_tables && data.rate_can_be_tabulated) {
        evaluate_tabulated_rate<rate>(state, rates);
    }
    else {
        evaluate_analytical_rate<rate>(state, rates);
    }

    // Set the density dependence

    apply_density_scaling<rate>(state, rates);

#ifdef SCREENING
    // Screen

    apply_screening<rate>(state, rates);
#endif

    // Branching ratios

    apply_branching<rate>(rates);
}

template<int rate, class Arr>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void fill_additional_rates (const Arr& intermediate_rates, rate_t& rates1, rate_t& rates2, rate_t& rates3)
{
    constexpr rhs_t data = RHS::rhs_data(rate);

    constexpr int rate1 = data.additional_reaction_1;
    constexpr int rate2 = data.additional_reaction_2;
    constexpr int rate3 = data.additional_reaction_3;

    if constexpr (rate1 >= 0) {
        constexpr int index = locate_intermediate_rate_index(rate1);

        static_assert(index >= 1 && index <= num_intermediate_reactions());

        rates1 = intermediate_rates(index);
    }

    if constexpr (rate2 >= 0) {
        constexpr int index = locate_intermediate_rate_index(rate2);

        static_assert(index >= 1 && index <= num_intermediate_reactions());

        rates2 = intermediate_rates(index);
    }

    if constexpr (rate3 >= 0) {
        constexpr int index = locate_intermediate_rate_index(rate3);

        static_assert(index >= 1 && index <= num_intermediate_reactions());

        rates3 = intermediate_rates(index);
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

// The components of ydot can either be the actual RHS terms (neqs)
// or separate tracking of positive and negative contributions to
// the RHS (2 * neqs, with the positive followed by the negative
// for each term in ydot).
template<int nrhs>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs (burn_t& burn_state, Array1D<Real, 1, nrhs>& ydot)
{
    static_assert(nrhs == neqs || nrhs == 2 * neqs);

    rhs_state_t rhs_state;

    rhs_state.rho = burn_state.rho;
    rhs_state.eta = burn_state.eta;
    rhs_state.y_e = burn_state.y_e;

    // Convert X to Y.
    for (int n = 1; n <= NumSpec; ++n) {
        rhs_state.y(n) = burn_state.xn[n-1] * aion_inv[n-1];
    }

    // Set up the state data, which is the same for all screening factors.
    fill_plasma_state(rhs_state.pstate, burn_state.T, burn_state.rho, rhs_state.y);

    // Initialize the rate temperature term.
    rhs_state.tf = get_tfactors(burn_state.T);
    if (use_tables) {
        rhs_state.tab.initialize(burn_state.T);
    }

    // Initialize the RHS terms.
    for (int n = 1; n <= nrhs; ++n) {
        ydot(n) = 0.0;
    }

    // Count up number of intermediate rates (rates that are used in any other reaction).
    constexpr int num_intermediate = num_intermediate_reactions();

    // We cannot have a zero-sized array, so just set the array size to 1 in that case.
    constexpr int intermediate_array_size = num_intermediate > 0 ? num_intermediate : 1;

    // Define forward and reverse (and d/dT) rate arrays.
    Array1D<rate_t, 1, intermediate_array_size> intermediate_rates;

    // Fill all intermediate rates first.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n)
    {
        constexpr int rate = n;

        constexpr int index = locate_intermediate_rate_index(rate);

        if constexpr (index >= 1) {
            construct_rate<rate>(rhs_state, intermediate_rates(index));
         }
    });

    // Loop over all rates, and then loop over all species, and if the
    // rate affects that given species, add its contribution to the RHS.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n1)
    {
        constexpr int rate = n1;

        rate_t rates;

        // We only need to compute the rate at this point if it's not intermediate. If it
        // is intermediate, retrieve it from the cached array.

        constexpr int index = locate_intermediate_rate_index(rate);
        if constexpr (index < 0) {
            construct_rate<rate>(rhs_state, rates);
        }
        else {
            rates = intermediate_rates(index);
        }

        // Locate all intermediate rates needed to augment this reaction.
        // To keep the problem bounded we assume that there are no more than
        // three intermediate reactions needed.

        rate_t rates1, rates2, rates3;

        fill_additional_rates<rate>(intermediate_rates, rates1, rates2, rates3);

        // Perform rate postprocessing, using additional reactions as inputs.
        // If there is no postprocessing for this rate, this will be a no-op.

        postprocess_rate<rate>(rhs_state, rates, rates1, rates2, rates3);

        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int species = n2;

            if constexpr (is_rate_used<species, rate>()) {
                constexpr int use_T_derivatives = 0;
                auto [forward_term, reverse_term] = rhs_term<species, rate, use_T_derivatives>(burn_state, rates);

                if constexpr (nrhs == 2 * neqs) {
                    if (forward_term >= 0.0_rt) {
                        ydot(2 * species - 1) += forward_term;
                        ydot(2 * species    ) += -reverse_term;
                    }
                    else {
                        ydot(2 * species - 1) += reverse_term;
                        ydot(2 * species    ) += -forward_term;
                    }
                }
                else {
                    ydot(species) += forward_term + reverse_term;
                }
            }
        });
    });

    // Evaluate the neutrino cooling.
#ifdef NEUTRINOS
    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(burn_state.T, burn_state.rho, burn_state.abar, burn_state.zbar,
           sneut, dsneutdt, dsneutdd, snuda, snudz);
#else
    Real sneut = 0.0;
#endif

    // Compute the energy RHS term.
    if constexpr (nrhs == 2 * neqs) {
        ydot(2 * net_ienuc - 1) = 0.0;
        ydot(2 * net_ienuc    ) = sneut;
    }
    else {
        ydot(net_ienuc) = -sneut;
    }

    constexpr_for<1, NumSpec+1>([&] (auto n)
    {
        constexpr int species = n;

        if constexpr (nrhs == 2 * neqs) {
            ydot(2 * net_ienuc - 1) += ener_gener_rate<species>(ydot(2 * species - 1));
            ydot(2 * net_ienuc    ) += ener_gener_rate<species>(ydot(2 * species    ));
        }
        else {
            ydot(net_ienuc) += ener_gener_rate<species>(ydot(species));
        }
    });
}

// Analytical Jacobian
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jac (burn_t& burn_state, ArrayUtil::MathArray2D<1, neqs, 1, neqs>& jac)
{
    rhs_state_t rhs_state;

    rhs_state.rho = burn_state.rho;
    rhs_state.eta = burn_state.eta;
    rhs_state.y_e = burn_state.y_e;

    // Convert X to Y.
    for (int n = 1; n <= NumSpec; ++n) {
        rhs_state.y(n) = burn_state.xn[n-1] * aion_inv[n-1];
    }

    // Set up the state data, which is the same for all screening factors.
    fill_plasma_state(rhs_state.pstate, burn_state.T, burn_state.rho, rhs_state.y);

    // Initialize the rate temperature term.
    rhs_state.tf = get_tfactors(burn_state.T);
    if (use_tables) {
        rhs_state.tab.initialize(burn_state.T);
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
    Array1D<rate_t, 1, intermediate_array_size> intermediate_rates;

    rate_t rates1, rates2, rates3;

    // Fill all intermediate rates first.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n)
    {
        constexpr int rate = n;

        constexpr int index = locate_intermediate_rate_index(rate);

        if constexpr (index >= 1) {
            construct_rate<rate>(rhs_state, intermediate_rates(index));
         }
    });

    // Loop over rates and compute Jacobian terms.
    constexpr_for<1, Rates::NumRates+1>([&] (auto n1)
    {
        constexpr int rate = n1;

        rate_t rates;

        // We only need to compute the rate at this point if it's not intermediate. If it
        // is intermediate, retrieve it from the cached array.

        constexpr int index = locate_intermediate_rate_index(rate);
        if constexpr (index < 0) {
            construct_rate<rate>(rhs_state, rates);
        }
        else {
            rates = intermediate_rates(index);
        }

        // Locate all intermediate rates needed to augment this reaction.
        // To keep the problem bounded we assume that there are no more than
        // three intermediate reactions needed.

        fill_additional_rates<rate>(intermediate_rates, rates1, rates2, rates3);

        // Perform rate postprocessing, using additional reactions as inputs.
        // If there is no postprocessing for this rate, this will be a no-op.

        postprocess_rate<rate>(rhs_state, rates, rates1, rates2, rates3);

        // Species Jacobian elements with respect to other species.
        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            [[maybe_unused]] constexpr int spec1 = n2;

            constexpr_for<1, NumSpec+1>([&] (auto n3)
            {
                [[maybe_unused]] constexpr int spec2 = n3;

                if constexpr (is_rate_used<spec1, rate>()) {
                    jac(spec1, spec2) += jac_term<spec1, spec2, rate>(burn_state, rates);
                }
            });
        });

        // Evaluate the Jacobian elements with respect to temperature.
        // We'll convert them from d/dT to d/de later.
        constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            [[maybe_unused]] constexpr int species = n2;

            if constexpr (is_rate_used<species, rate>()) {
                constexpr int use_T_derivatives = 1;
                auto [forward_term, reverse_term] = rhs_term<species, rate, use_T_derivatives>(burn_state, rates);
                jac(species, net_ienuc) += forward_term + reverse_term;
            }
        });
    });

    // Evaluate the neutrino cooling.
#ifdef NEUTRINOS
    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(burn_state.T, burn_state.rho, burn_state.abar, burn_state.zbar,
           sneut, dsneutdt, dsneutdd, snuda, snudz);
#else
    Real sneut = 0.0, dsneutdt = 0.0, dsneutdd = 0.0, snuda = 0.0, snudz = 0.0;
#endif

    jac(net_ienuc, net_ienuc) = -temperature_to_energy_jacobian(burn_state, dsneutdt);

    constexpr_for<1, NumSpec+1>([&] (auto j)
    {
        constexpr int species = j;

        // Energy generation rate Jacobian elements with respect to species.
        Real b1 = (-burn_state.abar * burn_state.abar * snuda + (NetworkProperties::zion(species) - burn_state.zbar) * burn_state.abar * snudz);
        jac(net_ienuc, species) = -b1;

        constexpr_for<1, NumSpec+1>([&] (auto i)
        {
            constexpr int s = i;

            jac(net_ienuc, species) += ener_gener_rate<s>(jac(s, species));
        });

        // Convert previously computed terms from d/dT to d/de.
        jac(species, net_ienuc) = temperature_to_energy_jacobian(burn_state, jac(species, net_ienuc));

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
