#ifndef rhs_utilities_H
#define rhs_utilities_H

#include <network_utilities.H>
#include <actual_network.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Array1D<Real, 1, NumSpec> construct_species_rhs (Array1D<Real, 1, NumSpec> const& y, rate_t const& rr)
{
    Array1D<Real, 1, NumSpec> dydt;

    for (int spec = 1; spec <= NumSpec; ++spec) {

        Array1D<Real, 1, Rates::NumRates> a = {0.0};

        // Loop through rates, adding the ones we care about
        // to the array that we will sum.

        int a_ctr = 0;

        for (int rate = 1; rate <= Rates::NumRates; ++rate) {
            rhs_t rhs_data = RHS::rhs_data(spec, rate);

            if (std::abs(rhs_data.prefactor) > 0.0) {
                ++a_ctr;

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

                a(a_ctr) = term;
            }
        }

        dydt(spec) = esum<Rates::NumRates>(a);
    }

    return dydt;
}

#endif
