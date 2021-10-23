#ifndef rhs_utilities_H
#define rhs_utilities_H

#include <network_utilities.H>
#include <actual_network.H>

// Forward declarations

namespace RHS
{

template<int species, int rate>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int is_rate_used ();

template<int species, int j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real rhs_term (const burn_t& state, rate_fr_t const& rr);

template<int spec1, int spec2, int j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real jac_term (const burn_t& state, rate_fr_t const& rr);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_fr_t const& rr);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_t const& rr);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array2D<Real, 1, NumSpec, 1, NumSpec> species_jac (const burn_t& state, rate_fr_t const& rr);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array2D<Real, 1, NumSpec, 1, NumSpec> species_jac (const burn_t& state, rate_t const& rr);

} // namespace RHS



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



namespace RHS
{

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
    constexpr rhs_t rhs_data = RHS::rhs_data(rate);

    static_assert(species >= 1 && species <= NumSpec);
    static_assert(rate >= 1 && rate <= Rates::NumRates);

    if (rhs_data.species_A == species ||
        rhs_data.species_B == species ||
        rhs_data.species_C == species ||
        rhs_data.species_D == species ||
        rhs_data.species_E == species ||
        rhs_data.species_F == species) {
        // Exclude intermediate rates that don't appear in the RHS.
        // We can identify these by the presence of an "extra" species
        // whose numerical ID is > NumSpec.
        if (rhs_data.species_A > NumSpec ||
            rhs_data.species_B > NumSpec ||
            rhs_data.species_C > NumSpec ||
            rhs_data.species_D > NumSpec ||
            rhs_data.species_E > NumSpec ||
            rhs_data.species_F > NumSpec) {
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
    constexpr rhs_t rhs_data = RHS::rhs_data(rate);

    int exponent = 0;

    if constexpr (rhs_data.species_A >= 0) {
        exponent += rhs_data.exponent_A;
    }
    if constexpr (rhs_data.species_B >= 0) {
        exponent += rhs_data.exponent_B;
    }
    if constexpr (rhs_data.species_C >= 0) {
        exponent += rhs_data.exponent_C;
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
    constexpr rhs_t rhs_data = RHS::rhs_data(rate);

    int exponent = 0;

    if constexpr (rhs_data.species_D >= 0) {
        exponent += rhs_data.exponent_D;
    }
    if constexpr (rhs_data.species_E >= 0) {
        exponent += rhs_data.exponent_E;
    }
    if constexpr (rhs_data.species_F >= 0) {
        exponent += rhs_data.exponent_F;
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
    constexpr int forward_exponent = RHS::density_exponent_forward<rate>();
    constexpr int reverse_exponent = RHS::density_exponent_reverse<rate>();

    Real density_term_forward = std::pow(rho, forward_exponent);
    Real density_term_reverse = std::pow(rho, reverse_exponent);

    fr *= density_term_forward;
    frdt *= density_term_forward;
    rr *= density_term_reverse;
    rrdt *= density_term_reverse;
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
constexpr Real rhs_term (const burn_t& state, rate_fr_t const& rr)
{
    constexpr rhs_t rhs_data = RHS::rhs_data(rate);

    // First, compute the Y * rate component of both the forward and
    // reverse reactions, which is the same regardless of which species
    // we're producing or consuming.

    Real forward_term = rr.forward(rate);

    if constexpr (rhs_data.species_A >= 0) {
        Real Y_A = state.xn[rhs_data.species_A-1] * aion_inv[rhs_data.species_A-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(rhs_data.exponent_A);
        forward_term *= std::pow(Y_A, rhs_data.exponent_A) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (rhs_data.species_B >= 0) {
        Real Y_B = state.xn[rhs_data.species_B-1] * aion_inv[rhs_data.species_B-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(rhs_data.exponent_B);
        forward_term *= std::pow(Y_B, rhs_data.exponent_B) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (rhs_data.species_C >= 0) {
        Real Y_C = state.xn[rhs_data.species_C-1] * aion_inv[rhs_data.species_C-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(rhs_data.exponent_C);
        forward_term *= std::pow(Y_C, rhs_data.exponent_C) /
                        static_cast<Real>(identical_particle_factor);
    }

    Real reverse_term = rr.reverse(rate);

    if constexpr (rhs_data.species_D >= 0) {
        Real Y_D = state.xn[rhs_data.species_D-1] * aion_inv[rhs_data.species_D-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(rhs_data.exponent_D);
        reverse_term *= std::pow(Y_D, rhs_data.exponent_D) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (rhs_data.species_E >= 0) {
        Real Y_E = state.xn[rhs_data.species_E-1] * aion_inv[rhs_data.species_E-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(rhs_data.exponent_E);
        reverse_term *= std::pow(Y_E, rhs_data.exponent_E) /
                        static_cast<Real>(identical_particle_factor);
    }

    if constexpr (rhs_data.species_F >= 0) {
        Real Y_F = state.xn[rhs_data.species_F-1] * aion_inv[rhs_data.species_F-1];
        constexpr int identical_particle_factor = RHS_impl::factorial(rhs_data.exponent_F);
        reverse_term *= std::pow(Y_F, rhs_data.exponent_F) /
                        static_cast<Real>(identical_particle_factor);
    }

    // Now compute the total contribution to this species.

    Real term = 0.0_rt;

    if constexpr (rhs_data.species_A == species) {
        term = rhs_data.number_A * (reverse_term - forward_term);
    }
    if constexpr (rhs_data.species_B == species) {
        term = rhs_data.number_B * (reverse_term - forward_term);
    }
    if constexpr (rhs_data.species_C == species) {
        term = rhs_data.number_C * (reverse_term - forward_term);
    }
    if constexpr (rhs_data.species_D == species) {
        term = rhs_data.number_D * (forward_term - reverse_term);
    }
    if constexpr (rhs_data.species_E == species) {
        term = rhs_data.number_E * (forward_term - reverse_term);
    }
    if constexpr (rhs_data.species_F == species) {
        term = rhs_data.number_F * (forward_term - reverse_term);
    }

    return term;
}

// Calculate the array of RHS terms over all species.
// Legacy interface using old-style rate array.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_t const& rr)
{
    rate_fr_t rr_fr;

    // Convert old-style rates to new-style
    for (int j = 1; j <= Rates::NumRatesFR; ++j) {
        rr_fr.forward(j) = rr.rates(2 * j - 1);
        rr_fr.reverse(j) = rr.rates(2 * j);
    }

    return species_rhs(state, rr_fr);
}

// Calculate the array of RHS terms over all species.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_fr_t const& rr)
{
    Array1D<Real, 1, NumSpec> spec_rhs = {0.0};

    // For every species, loop over all rates, and if that rate involves
    // the given species, add its contribution to the RHS.

    RHS::constexpr_for<1, NumSpec+1>([&] (auto n1)
    {
        constexpr int species = n1;
        RHS::constexpr_for<1, Rates::NumRatesFR+1>([&] (auto n2)
        {
            constexpr int rate = n2;
            if constexpr (is_rate_used<species, rate>()) {
                spec_rhs(species) += RHS::rhs_term<species, rate>(state, rr);
            }
        });
    });

    return spec_rhs;
}

// Calculate the energy generation rate, including neutrino losses.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real ener_rhs (const burn_t& state, Array1D<Real, 1, NumSpec>& dydt)
{
    using namespace network;

    Real Xdot = 0.0_rt;

    // Sum the mass fraction time derivatives

    for (int i = 1; i <= NumSpec; ++i) {
        Xdot += dydt(i) * mion(i);
    }

    // This is basically e = m c**2

    Real dedt = Xdot * C::Legacy::enuc_conv2;

    // Get the neutrino losses

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    dedt -= sneut;

    return dedt;
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
constexpr Real jac_term (const burn_t& state, rate_fr_t const& rr)
{
    constexpr rhs_t rhs_data = RHS::rhs_data(rate);

    Real forward_term = 0.0_rt;

    if constexpr (spec2 == rhs_data.species_A || spec2 == rhs_data.species_B || spec2 == rhs_data.species_C) {

        forward_term = rr.forward(rate);

        if constexpr (rhs_data.species_A >= 0) {
            Real Y_A = state.xn[rhs_data.species_A-1] * aion_inv[rhs_data.species_A-1];

            constexpr int exponent = rhs_data.exponent_A;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            // Compute the forward term. It only has a contribution if
            // the species we're taking the derivative with respect to
            // (spec2) is one of (A, B, C). For the species that is spec2,
            // we take the derivative by multiplying by the current exponent
            // and then decrementing the exponent in the term.

            int exp = exponent;

            if constexpr (spec2 == rhs_data.species_A) {
                forward_term *= exp;
                --exp;
            }

            forward_term *= std::pow(Y_A, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (rhs_data.species_B >= 0) {
            Real Y_B = state.xn[rhs_data.species_B-1] * aion_inv[rhs_data.species_B-1];

            constexpr int exponent = rhs_data.exponent_B;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == rhs_data.species_B) {
                forward_term *= exp;
                --exp;
            }

            forward_term *= std::pow(Y_B, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (rhs_data.species_C >= 0) {
            Real Y_C = state.xn[rhs_data.species_C-1] * aion_inv[rhs_data.species_C-1];

            constexpr int exponent = rhs_data.exponent_C;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == rhs_data.species_C) {
                forward_term *= exp;
                --exp;
            }

            forward_term *= std::pow(Y_C, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

    }

    Real reverse_term = 0.0_rt;

    if constexpr (spec2 == rhs_data.species_D || spec2 == rhs_data.species_E || spec2 == rhs_data.species_F) {

        reverse_term = rr.reverse(rate);

        if constexpr (rhs_data.species_D >= 0) {
            Real Y_D = state.xn[rhs_data.species_D-1] * aion_inv[rhs_data.species_D-1];

            constexpr int exponent = rhs_data.exponent_D;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == rhs_data.species_D) {
                reverse_term *= exp;
                --exp;
            }

            reverse_term *= std::pow(Y_D, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (rhs_data.species_E >= 0) {
            Real Y_E = state.xn[rhs_data.species_E-1] * aion_inv[rhs_data.species_E-1];

            constexpr int exponent = rhs_data.exponent_E;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == rhs_data.species_E) {
                reverse_term *= exp;
                --exp;
            }

            reverse_term *= std::pow(Y_E, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

        if constexpr (rhs_data.species_F >= 0) {
            Real Y_F = state.xn[rhs_data.species_F-1] * aion_inv[rhs_data.species_F-1];

            constexpr int exponent = rhs_data.exponent_F;

            constexpr int identical_particle_factor = RHS_impl::factorial(exponent);

            int exp = exponent;

            if constexpr (spec2 == rhs_data.species_F) {
                reverse_term *= exp;
                --exp;
            }

            reverse_term *= std::pow(Y_F, exp) /
                            static_cast<Real>(identical_particle_factor);
        }

    }

    // Now compute the total contribution to this species.

    Real term = 0.0_rt;

    if constexpr (rhs_data.species_A == spec1) {
        term = rhs_data.number_A * (reverse_term - forward_term);
    }
    if constexpr (rhs_data.species_B == spec1) {
        term = rhs_data.number_B * (reverse_term - forward_term);
    }
    if constexpr (rhs_data.species_C == spec1) {
        term = rhs_data.number_C * (reverse_term - forward_term);
    }
    if constexpr (rhs_data.species_D == spec1) {
        term = rhs_data.number_D * (forward_term - reverse_term);
    }
    if constexpr (rhs_data.species_E == spec1) {
        term = rhs_data.number_E * (forward_term - reverse_term);
    }
    if constexpr (rhs_data.species_F == spec1) {
        term = rhs_data.number_F * (forward_term - reverse_term);
    }

    return term;
}

// Calculate the Jacobian terms that are derivatives of one species
// with respect to another.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array2D<Real, 1, NumSpec, 1, NumSpec> species_jac (const burn_t& state, rate_t const& rr)
{
    rate_fr_t rr_fr;

    // Convert old-style rates to new-style
    for (int j = 1; j <= Rates::NumRatesFR; ++j) {
        rr_fr.forward(j) = rr.rates(2 * j - 1);
        rr_fr.reverse(j) = rr.rates(2 * j);
    }

    return RHS::species_jac(state, rr_fr);
}

// Calculate the Jacobian terms that are derivatives of one species
// with respect to another.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array2D<Real, 1, NumSpec, 1, NumSpec> species_jac (const burn_t& state, rate_fr_t const& rr)
{
    Array2D<Real, 1, NumSpec, 1, NumSpec> spec_jac = {0.0};

    // For every spec1 and spec2 pair, loop over all rates, and if that rate
    // involves spec1, add its contribution to the Jacobian.

    RHS::constexpr_for<1, NumSpec+1>([&] (auto n1)
    {
        constexpr int spec1 = n1;
        RHS::constexpr_for<1, NumSpec+1>([&] (auto n2)
        {
            constexpr int spec2 = n2;
            RHS::constexpr_for<1, Rates::NumRatesFR+1>([&] (auto n3)
            {
                constexpr int rate = n3;
                if constexpr (is_rate_used<spec1, rate>()) {
                    spec_jac(spec1, spec2) += RHS::jac_term<spec1, spec2, rate>(state, rr);
                }
            });
        });
    });

    return spec_jac;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, 6> equilibrium_state (const rhs_t& rhs_data, const burn_t& state,
                                       Real forward_rate, Real reverse_rate)
{
    int n_A = rhs_data.number_A;
    int n_B = rhs_data.number_B;
    int n_C = rhs_data.number_C;
    int n_D = rhs_data.number_D;
    int n_E = rhs_data.number_E;
    int n_F = rhs_data.number_F;

    int e_A = rhs_data.exponent_A;
    int e_B = rhs_data.exponent_B;
    int e_C = rhs_data.exponent_C;
    int e_D = rhs_data.exponent_D;
    int e_E = rhs_data.exponent_E;
    int e_F = rhs_data.exponent_F;

    Real Y_A = 0.0_rt;
    if (rhs_data.species_A >= 0) {
        Y_A = state.xn[rhs_data.species_A-1] * aion_inv[rhs_data.species_A-1];
    }

    Real Y_B = 0.0_rt;
    if (rhs_data.species_B >= 0) {
        Y_B = state.xn[rhs_data.species_B-1] * aion_inv[rhs_data.species_B-1];
    }

    Real Y_C = 0.0_rt;
    if (rhs_data.species_C >= 0) {
        Y_C = state.xn[rhs_data.species_C-1] * aion_inv[rhs_data.species_C-1];
    }

    Real Y_D = 0.0_rt;
    if (rhs_data.species_D >= 0) {
        Y_D = state.xn[rhs_data.species_D-1] * aion_inv[rhs_data.species_D-1];
    }

    Real Y_E = 0.0_rt;
    if (rhs_data.species_E >= 0) {
        Y_E = state.xn[rhs_data.species_E-1] * aion_inv[rhs_data.species_E-1];
    }

    Real Y_F = 0.0_rt;
    if (rhs_data.species_F >= 0) {
        Y_F = state.xn[rhs_data.species_F-1] * aion_inv[rhs_data.species_F-1];
    }

    Real k_f = forward_rate;
    Real k_r = reverse_rate;

    Real Y_A_eq = 0.0_rt;
    Real Y_B_eq = 0.0_rt;
    Real Y_C_eq = 0.0_rt;
    Real Y_D_eq = 0.0_rt;
    Real Y_E_eq = 0.0_rt;
    Real Y_F_eq = 0.0_rt;

    // The general form of our reactions has forward reactions of the form:
    // N * (Y_A**e_A * Y_B**e_B * Y_C**e_C * k_f)
    // and reverse reactions of the form:
    // N * (Y_D**e_D * Y_B**e_E * Y_F**e_F * k_r)
    // where N is the number of particles consumed/produced, Y is the molar
    // fraction, and e is an integer exponent coefficient, k_f is the forward
    // rate, and k_r is the reverse rate.
    //
    // If the reaction is in equilibrium, we will satisfy
    // Y_A**e_A * Y_B**e_B * Y_C**e_C * k_f = Y_D**e_D * Y_B**e_E * Y_F**e_F * k_r
    // We expect that often this relation will not be satisfied exactly,
    // but sufficiently close that we will likely achieve equilibrium soon.
    //
    // So, given a state and a set of rates, we want to determine the corresponding
    // equilibrium solution. This involves solving the t -> infinity limit of the ODE.
    // For certain common values of the exponents this solution is easy to write down
    // analytically. We will solve for these cases specifically and ignore other cases.

    if (e_A == 1 && e_B == 1 && e_C == 0 && e_D == 1 && e_E == 0 && e_F == 0 &&
        n_A == 1 && n_B == 1 && n_C == 0 && n_D == 1 && n_E == 0 && n_F == 0)
    {
        // This reaction looks like
        // A + B --> D
        // and the corresponding ODE for species A is
        // dY_A / dt = -k_f * Y_A * Y_B + k_r * Y_D
        //
        // The conserved scalars are
        // c_1 = Y_A - Y_B
        // c_2 = Y_A + Y_D
        //
        // The species A ODE in terms of only Y_A is
        // dY_A / dt = a * Y_A**2 + b * Y_A + c
        // where
        // a = -k_f
        // b = c_1 * k_f - k_r
        // c = k_r * c_2
        //
        // The equilibrium solution for species A is
        // Y_A_eq = (-1 / 2a) * (b + sqrt(b * b - 4 * a * c))
        //
        // The equilibrium solutions for the other species are
        // Y_B_eq = Y_A_eq - c_1
        // Y_D_eq = c_2 - Y_A_eq

        Real c_1 = Y_A - Y_B;
        Real c_2 = Y_A + Y_D;

        Real a = -k_f;
        Real b = c_1 * k_f - k_r;
        Real c = k_r * c_2;

        Y_A_eq = (-0.5_rt / a) * (b + std::sqrt(b * b - 4.0_rt * a * c));
        Y_B_eq = Y_A_eq - c_1;
        Y_D_eq = c_2 - Y_A_eq;
    }

    Real X_A_eq = 0.0_rt;
    if (rhs_data.species_A >= 0) {
        X_A_eq = Y_A_eq * aion[rhs_data.species_A-1];
    }

    Real X_B_eq = 0.0_rt;
    if (rhs_data.species_B >= 0) {
        X_B_eq = Y_B_eq * aion[rhs_data.species_B-1];
    }

    Real X_C_eq = 0.0_rt;
    if (rhs_data.species_C >= 0) {
        X_C_eq = Y_C_eq * aion[rhs_data.species_C-1];
    }

    Real X_D_eq = 0.0_rt;
    if (rhs_data.species_D >= 0) {
        X_D_eq = Y_D_eq * aion[rhs_data.species_D-1];
    }

    Real X_E_eq = 0.0_rt;
    if (rhs_data.species_E >= 0) {
        X_E_eq = Y_E_eq * aion[rhs_data.species_E-1];
    }

    Real X_F_eq = 0.0_rt;
    if (rhs_data.species_F >= 0) {
        X_F_eq = Y_F_eq * aion[rhs_data.species_F-1];
    }

    Array1D<Real, 1, 6> eq_state = {X_A_eq, X_B_eq, X_C_eq, X_D_eq, X_E_eq, X_F_eq};

    return eq_state;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
bool is_reaction_in_equilibrium (const rhs_t& rhs_data, const burn_t& state,
                                 const Array1D<Real, 1, 6>& equilibrium_state)
{
    // Check whether all species are within a given certain tolerance of the
    // corresponding equilibrium mass fractions.

    const Real equilibrium_tolerance = 0.005_rt;

    bool equilibrium = true;

    if (rhs_data.species_A >= 0) {
        Real X = state.xn[rhs_data.species_A - 1];
        Real X_eq = equilibrium_state(1);

        if (std::abs(X - X_eq) > equilibrium_tolerance * X || X < atol_spec) {
            equilibrium = false;
        }
    }

    if (rhs_data.species_B >= 0) {
        Real X = state.xn[rhs_data.species_B - 1];
        Real X_eq = equilibrium_state(2);

        if (std::abs(X - X_eq) > equilibrium_tolerance * X || X < atol_spec) {
            equilibrium = false;
        }
    }

    if (rhs_data.species_C >= 0) {
        Real X = state.xn[rhs_data.species_C - 1];
        Real X_eq = equilibrium_state(3);

        if (std::abs(X - X_eq) > equilibrium_tolerance * X || X < atol_spec) {
            equilibrium = false;
        }
    }

    if (rhs_data.species_D >= 0) {
        Real X = state.xn[rhs_data.species_D - 1];
        Real X_eq = equilibrium_state(4);

        if (std::abs(X - X_eq) > equilibrium_tolerance * X || X < atol_spec) {
            equilibrium = false;
        }
    }

    if (rhs_data.species_E >= 0) {
        Real X = state.xn[rhs_data.species_E - 1];
        Real X_eq = equilibrium_state(5);

        if (std::abs(X - X_eq) > equilibrium_tolerance * X || X < atol_spec) {
            equilibrium = false;
        }
    }

    if (rhs_data.species_F >= 0) {
        Real X = state.xn[rhs_data.species_F - 1];
        Real X_eq = equilibrium_state(6);

        if (std::abs(X - X_eq) > equilibrium_tolerance * X || X < atol_spec) {
            equilibrium = false;
        }
    }

    return equilibrium;
}

} // namespace RHS

#endif
