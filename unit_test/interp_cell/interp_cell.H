#ifndef INTERP_CELL_H
#define INTERP_CELL_H

#include <AMReX_Array.H>

#include <extern_parameters.H>

#include <interp_tools.H>

amrex::Real test_function(amrex::Real x) {
    amrex::Real x0{-2.5};
    return (2.0_rt/3.0_rt) * amrex::Math::powi<3>(x - x0) - 1.5_rt * amrex::Math::powi<2>(x - x0) + 7.0_rt * (x - x0) - 10.0_rt;
}

AMREX_INLINE
void interp_cell_c()
{

    // create some fake data
    // our test function is a cubic so we should fit it
    // exactly with out interpolation
    amrex::Array1D<amrex::Real, 1, 10> x_array{-3.2, -2.5, -1.23, -0.6, 0.01, 0.54, 0.98, 1.35, 1.87, 2.5};
    amrex::Array1D<amrex::Real, 1, 10> f_array{};
    std::cout << "initial data: " << std::endl;
    for (int i = x_array.lo(); i <= x_array.hi(); ++i) {
        f_array(i) = test_function(x_array(i));
        std::cout << x_array(i) << ", " << f_array(i) << std::endl;
    }
    std::cout << std::endl;

    // now we see if we recover things well
    {
        amrex::Real x0{-1.1};
        constexpr int do_derivative = 1;
        const auto [f, fprime] = interp_net::cubic_interp_uneven<do_derivative>(x0, x_array, f_array);

        std::cout << f << " " << test_function(x0) << std::endl;
    }

    {
        amrex::Real x0{1.42};
        constexpr int do_derivative = 1;
        const auto [f, fprime] = interp_net::cubic_interp_uneven<do_derivative>(x0, x_array, f_array);

        std::cout << f << " " << test_function(x0) << std::endl;
    }




}
#endif
