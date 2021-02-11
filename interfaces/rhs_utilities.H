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
constexpr Real rhs_term (const burn_t& state, rate_t const& rr);

template<int species, int N>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, N>
make_RHS_Array1D (const burn_t& state, rate_t const& rr);

template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real species_rhs_n (const burn_t& state, rate_t const& rr);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_t const& rr);

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
make_RHS_Array1D_impl (const burn_t& state, rate_t const& rr, std::integer_sequence<int, j...>)
{
    return {{RHS::rhs_term<species, j+1>(state, rr)...}};
}

// Implicitly construct an Array1D by expanding the integer sequence.
// Note that the integer sequence is zero-indexed but the terms are
// one-indexed.
template<int... species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs_impl (const burn_t& state, rate_t const& rr, std::integer_sequence<int, species...>)
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
    constexpr rhs_t rhs_data = RHS::rhs_data(species, rate);

    static_assert(species >= 1 && species <= NumSpec);
    static_assert(rate >= 1 && rate <= Rates::NumRates);

    if (std::abs(rhs_data.prefactor) > 0.0) {
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
// It is of the form K * Y(1) * Y(2) * Y(3) * rate,
// where K is a prefactor constant, rate is the
// reaction rate, and Y(1), Y(2), and Y(3) are up to
// three molar fractions that participate (for one-body,
// two-body, and three-body reactions, respectively). If
// a given reaction uses fewer than three bodies, we infer
// this by calling its index -1 and then not accessing it
// in the multiplication.
template<int species, int j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real rhs_term (const burn_t& state, rate_t const& rr)
{
    constexpr int rate = RHS_impl::rhs_term_impl<species, j>(std::make_integer_sequence<int, Rates::NumRates>{});

    constexpr rhs_t rhs_data = RHS::rhs_data(species, rate);

    Real term = rhs_data.prefactor;

    if (rhs_data.specindex1 >= 0) {
        term *= state.xn[rhs_data.specindex1-1] * aion_inv[rhs_data.specindex1-1];
    }

    if (rhs_data.specindex2 >= 0) {
        term *= state.xn[rhs_data.specindex2-1] * aion_inv[rhs_data.specindex2-1];
    }

    if (rhs_data.specindex3 >= 0) {
        term *= state.xn[rhs_data.specindex3-1] * aion_inv[rhs_data.specindex3-1];
    }

    term *= rr.rates(rate);

    return term;
}

// Calculate the set of RHS terms.
template<int species, int N>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, N>
make_RHS_Array1D (const burn_t& state, rate_t const& rr)
{
    return RHS_impl::make_RHS_Array1D_impl<species>(state, rr, std::make_integer_sequence<int, N>{});
}

// Calculate the RHS for a given species by constructing the array of terms
// and then summing them up.
template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real species_rhs_n (const burn_t& state, rate_t const& rr)
{
    constexpr int nrhs = num_rhs<species>();
    Array1D<Real, 1, nrhs> a = make_RHS_Array1D<species, nrhs>(state, rr);
    return esum<nrhs>(a);
}

// Calculate the array of RHS terms over all species.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_t const& rr)
{
    return RHS_impl::species_rhs_impl(state, rr, std::make_integer_sequence<int, NumSpec>{});
}

} // namespace RHS

#endif
