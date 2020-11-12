#ifndef rhs_utilities_H
#define rhs_utilities_H

#include <actual_network.H>
#include <network_utilities.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr int num_rhs (int species)
{
    int nrhs = 0;

    // Query how many terms are non-zero. At most
    // the number of non-zero terms would be the
    // number of rates, so we have an upper limit
    // to loop over.
    for (int n = 0; n < Rates::NumRates; ++n) {
        rhs_t rhs_data = RHS::rhs_data(species, n);
        if (std::abs(rhs_data.prefactor) > 0.0) {
            ++nrhs;
        }
    }

    return nrhs;
}

template<int spec, int j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Real rhs_term (Array1D<Real, 1, NumSpec> const& y, rate_t const& rr)
{
    // Obtain the j'th term in the RHS (one-indexed).
    // We loop through the rates in order and
    // count up the order in which the rates are
    // used in the RHS (determined by the prefactor
    // being non-zero). We then return the rate
    // index of that term.

    rhs_t rhs_data{};

    int rate;
    int num_rates_used = 0;
    for (rate = 0; rate < Rates::NumRates; ++rate) {
        rhs_data = RHS::rhs_data(spec, rate);
        if (std::abs(rhs_data.prefactor) > 0.0) {
            ++num_rates_used;

            if (num_rates_used == j) {
                break;
            }
        }
    }

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

template<int spec, int... j>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, sizeof...(j)>
make_RHS_Array1D_impl(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr, std::integer_sequence<int, j...>) {
    return {{rhs_term<spec, j+1>(y, rr)...}};
}

template<int spec, int N>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
constexpr Array1D<Real, 1, N>
make_RHS_Array1D(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr) {
    return make_RHS_Array1D_impl<spec>(y, rr, std::make_integer_sequence<int, N>{});
}

template<int spec>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real spec_rhs(Array1D<Real, 1, NumSpec> const& y, rate_t const& rr) {
    constexpr int nrhs = num_rhs(spec);
    Array1D<Real, 1, nrhs> a = make_RHS_Array1D<spec, nrhs>(y, rr);
    return esum<nrhs>(a);
}

#endif
