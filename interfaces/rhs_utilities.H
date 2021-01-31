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

template<int species, int N>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, N>
make_RHS_Array1D (const burn_t& state, rate_fr_t const& rr);

template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real species_rhs_n (const burn_t& state, rate_fr_t const& rr);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_fr_t const& rr);

} // namespace RHS



// Internal implementation functions

namespace RHS_impl {

// Base case for recursion.
template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int add_nonzero_rates ()
{
    return 0;
}

// Recursively add nonzero rates.
template<int species, int rate, int... rates>
constexpr int add_nonzero_rates ()
{
    return RHS::is_rate_used<species, rate + 1>() + add_nonzero_rates<species, rates...>();
}

template<int species, int... rates>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int num_rhs_impl (std::integer_sequence<int, rates...>)
{
    return add_nonzero_rates<species, rates...>();
}

// Base case for recursion (no rates were found).
template<int species, int j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int actual_rhs_term_impl (int count)
{
    return -1;
}

// If we have counted up j nonzero rates, return the
// current rate index. Otherwise, keep counting.
template<int species, int j, int rate, int... rates>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int actual_rhs_term_impl (int count)
{
    count += RHS::is_rate_used<species, rate + 1>();

    if (count == j) {
        return rate + 1;
    }
    else {
        return actual_rhs_term_impl<species, j, rates...>(count);
    }
}

// Return the rate corresponding to the j'th contribution to the RHS.
// We obtain this by looping through the rates, and counting each one
// with a nonzero contribution. We stop when we have hit j of them.
template<int species, int j, int... rates>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int rhs_term_impl (std::integer_sequence<int, rates...>)
{
    return actual_rhs_term_impl<species, j, rates...>(0);
}

// Implicitly construct an Array1D by expanding the integer sequence.
// Note that the integer sequence is zero-indexed but the terms are
// one-indexed.
template<int species, int... j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, sizeof...(j)>
make_RHS_Array1D_impl (const burn_t& state, rate_fr_t const& rr, std::integer_sequence<int, j...>)
{
    return {{RHS::rhs_term<species, j+1>(state, rr)...}};
}

// Implicitly construct an Array1D by expanding the integer sequence.
// Note that the integer sequence is zero-indexed but the terms are
// one-indexed.
template<int... species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs_impl (const burn_t& state, rate_fr_t const& rr, std::integer_sequence<int, species...>)
{
    return {{RHS::species_rhs_n<species + 1>(state, rr)...}};
}

} // namespace RHS_impl



namespace RHS
{

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
        return 1;
    } else {
        return 0;
    }
}

// Count number of nonzero RHS terms.
template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int num_rhs ()
{
    return RHS_impl::num_rhs_impl<species>(std::make_integer_sequence<int, Rates::NumRates>{});
}

// Calculate the j'th RHS term for a given species.
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
// If a given reaction uses fewer than three species, we infer
// this by calling its index -1 and then not accessing it
// in the multiplication.
template<int species, int j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real rhs_term (const burn_t& state, rate_fr_t const& rr)
{
    constexpr int rate = RHS_impl::rhs_term_impl<species, j>(std::make_integer_sequence<int, Rates::NumRatesFR>{});

    constexpr rhs_t rhs_data = RHS::rhs_data(rate);

    // First, compute the Y * rate component of both the forward and
    // reverse reactions, which is the same regardless of which species
    // we're producing or consuming.

    Real forward_term = rr.forward_rates(rate);

    if constexpr (rhs_data.species_A >= 0) {
        Real Y_A = state.xn[rhs_data.species_A-1] * aion_inv[rhs_data.species_A-1];
        forward_term *= std::pow(Y_A, rhs_data.exponent_A_forward);
    }

    if constexpr (rhs_data.species_B >= 0) {
        Real Y_B = state.xn[rhs_data.species_B-1] * aion_inv[rhs_data.species_B-1];
        forward_term *= std::pow(Y_B, rhs_data.exponent_B_forward);
    }

    if constexpr (rhs_data.species_C >= 0) {
        Real Y_C = state.xn[rhs_data.species_C-1] * aion_inv[rhs_data.species_C-1];
        forward_term *= std::pow(Y_C, rhs_data.exponent_C_forward);
    }

    Real reverse_term = rr.reverse_rates(rate);

    if (rhs_data.species_D >= 0) {
        Real Y_D = state.xn[rhs_data.species_D-1] * aion_inv[rhs_data.species_D-1];
        reverse_term *= std::pow(Y_D, rhs_data.exponent_D_reverse);
    }

    if (rhs_data.species_E >= 0) {
        Real Y_E = state.xn[rhs_data.species_E-1] * aion_inv[rhs_data.species_E-1];
        reverse_term *= std::pow(Y_E, rhs_data.exponent_E_reverse);
    }

    if (rhs_data.species_F >= 0) {
        Real Y_F = state.xn[rhs_data.species_F-1] * aion_inv[rhs_data.species_F-1];
        reverse_term *= std::pow(Y_F, rhs_data.exponent_F_reverse);
    }

    // Now compute the total contribution to this species.

    Real term = 0.0_rt;

    if constexpr (rhs_data.species_A == species) {
        term = rhs_data.produced_A_reverse * reverse_term - rhs_data.consumed_A_forward * forward_term;
    }
    if constexpr (rhs_data.species_B == species) {
        term = rhs_data.produced_B_reverse * reverse_term - rhs_data.consumed_B_forward * forward_term;
    }
    if constexpr (rhs_data.species_C == species) {
        term = rhs_data.produced_C_reverse * reverse_term - rhs_data.consumed_C_forward * forward_term;
    }
    if constexpr (rhs_data.species_D == species) {
        term = rhs_data.produced_D_forward * forward_term - rhs_data.consumed_D_reverse * reverse_term;
    }
    if constexpr (rhs_data.species_E == species) {
        term = rhs_data.produced_E_forward * forward_term - rhs_data.consumed_E_reverse * reverse_term;
    }
    if constexpr (rhs_data.species_F == species) {
        term = rhs_data.produced_F_forward * forward_term - rhs_data.consumed_F_reverse * reverse_term;
    }

    return term;
}

// Calculate the set of RHS terms.
template<int species, int N>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, N>
make_RHS_Array1D (const burn_t& state, rate_fr_t const& rr)
{
    return RHS_impl::make_RHS_Array1D_impl<species>(state, rr, std::make_integer_sequence<int, N>{});
}

// Calculate the RHS for a given species by constructing the array of terms
// and then summing them up.
template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real species_rhs_n (const burn_t& state, rate_fr_t const& rr)
{
    constexpr int nrhs = num_rhs<species>();
    Array1D<Real, 1, nrhs> a = make_RHS_Array1D<species, nrhs>(state, rr);
    return esum<nrhs>(a);
}

// Calculate the array of RHS terms over all species.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_fr_t const& rr)
{
    return RHS_impl::species_rhs_impl(state, rr, std::make_integer_sequence<int, NumSpec>{});
}

// Overload of the above for the legacy rates.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_t const& rr)
{
    rate_fr_t rr_fr;

    // Convert old-style rates to new-style
    for (int j = 1; j <= Rates::NumRatesFR; ++j) {
        rr_fr.forward_rates(j) = rr.rates(2 * j - 1);
        rr_fr.reverse_rates(j) = rr.rates(2 * j);
    }

    return RHS_impl::species_rhs_impl(state, rr_fr, std::make_integer_sequence<int, NumSpec>{});
}

} // namespace RHS

#endif
