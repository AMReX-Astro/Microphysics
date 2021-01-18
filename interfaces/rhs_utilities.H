#ifndef rhs_utilities_H
#define rhs_utilities_H

#include <network_utilities.H>
#include <actual_network.H>

template<int species, int j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real rhs_term (const burn_t& state, rate_t const& rr);

template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real species_rhs_n (const burn_t& state, rate_t const& rr);

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_t const& rr);

namespace rhs_impl {

// Implicitly construct an Array1D by expanding the integer sequence.
// Note that the integer sequence is zero-indexed but the terms are
// one-indexed.
template<int species, int... j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, sizeof...(j)>
make_RHS_Array1D_impl (const burn_t& state, rate_t const& rr, std::integer_sequence<int, j...>)
{
    return {{rhs_term<species, j+1>(state, rr)...}};
}

// Calculate the set of RHS terms.
template<int species, int N>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, N>
make_RHS_Array1D (const burn_t& state, rate_t const& rr)
{
    return make_RHS_Array1D_impl<species>(state, rr, std::make_integer_sequence<int, N>{});
}

// Implicitly construct an Array1D by expanding the integer sequence.
// Note that the integer sequence is zero-indexed but the terms are
// one-indexed.
template<int... species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs_impl (const burn_t& state, rate_t const& rr, std::integer_sequence<int, species...>)
{
    return {{species_rhs_n<species + 1>(state, rr)...}};
}

} // namespace rhs_impl

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
    constexpr int rate = j;

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

// Calculate the RHS for a given species by constructing the array of terms
// and then summing them up.
template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real species_rhs_n (const burn_t& state, rate_t const& rr)
{
    constexpr int nrhs = Rates::NumRates;
    Array1D<Real, 1, nrhs> a = rhs_impl::make_RHS_Array1D<species, nrhs>(state, rr);
    return esum<nrhs>(a);
}

// Calculate the array of RHS terms over all species.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> species_rhs (const burn_t& state, rate_t const& rr)
{
    return rhs_impl::species_rhs_impl(state, rr, std::make_integer_sequence<int, NumSpec>{});
}

#endif
