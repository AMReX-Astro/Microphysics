#ifndef rhs_utilities_H
#define rhs_utilities_H

#include <actual_network.H>
#include <network_utilities.H>

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
    return is_rate_used<species, rate + 1>() + add_nonzero_rates<species, rates...>();
}

template<int species, int... rates>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int num_rhs_impl (std::integer_sequence<int, rates...>)
{
    return add_nonzero_rates<species, rates...>();
}

// Count number of nonzero RHS terms.
template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int num_rhs ()
{
    return num_rhs_impl<species>(std::make_integer_sequence<int, Rates::NumRates>{});
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
    count += is_rate_used<species, rate + 1>();

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
constexpr Real rhs_term (Array1D<Real, 1, NumSpec> const& y, rate_t const& rr)
{
    constexpr int rate = rhs_term_impl<species, j>(std::make_integer_sequence<int, Rates::NumRates>{});

    constexpr rhs_t rhs_data = RHS::rhs_data(species, rate);

    Real term = rhs_data.prefactor;

    if (rhs_data.specindex1 >= 0) {
        term *= y(rhs_data.specindex1);
    }

    if (rhs_data.specindex2 >= 0) {
        term *= y(rhs_data.specindex2);
    }

    if (rhs_data.specindex3 >= 0) {
        term *= y(rhs_data.specindex3);
    }

    term *= rr.rates(rate);

    return term;
}

template<int species, int... j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, sizeof...(j)>
make_RHS_Array1D_impl(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr, std::integer_sequence<int, j...>) {
    return {{rhs_term<species, j+1>(y, rr)...}};
}

template<int species, int N>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, N>
make_RHS_Array1D(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr) {
    return make_RHS_Array1D_impl<species>(y, rr, std::make_integer_sequence<int, N>{});
}

template<int species>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real species_rhs(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr) {
    constexpr int nrhs = num_rhs<species>();
    Array1D<Real, 1, nrhs> a = make_RHS_Array1D<species, nrhs>(y, rr);
    return esum<nrhs>(a);
}

#endif
