#include <catch2/catch.hpp>
#include "util/quadrature.H"
#include <AMReX_REAL.H>

TEST_CASE("Gauss-Kronrod integrates a polynomial exactly", "[quadrature]")
{
    // f(x) = x^2 + 2x + 3 in [a,b]
    auto f = [] AMREX_GPU_HOST_DEVICE (amrex::Real x) {
        return x*x + 2.0_rt*x + 3.0_rt;
    };
    amrex::Real a = -1.0_rt;
    amrex::Real b = 2.0_rt;
    // Analytical integral: x^3/3 + x^2 + 3x
    auto F = [&](amrex::Real x) {
        return x*x*x/3.0_rt + x*x + 3.0_rt*x;
    };
    amrex::Real exact = F(b) - F(a);
    auto [res, err] = gauss_kronrod_15_7(f, a, b);
    REQUIRE(err == Approx(0.0_rt).margin(1e-12));
    REQUIRE(res == Approx(exact).epsilon(1e-12));
}