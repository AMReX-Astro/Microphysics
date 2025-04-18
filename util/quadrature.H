#ifndef UTIL_QUADRATURE_H
#define UTIL_QUADRATURE_H

#include <AMReX_REAL.H>
#include <AMReX_Math.H>

/// \file quadrature.H
/// \brief Host/device inline Gauss–Kronrod (7/15) quadrature

/// Computes \f$\int_a^b f(x)\,dx\f$ using the 15-point Kronrod rule
/// with embedded 7-point Gauss estimate. Returns {integral, error}.
/// F must be a GPU-friendly functor: Real operator()(Real x) const;
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::pair<amrex::Real,amrex::Real>
gauss_kronrod_15_7(F f, amrex::Real a, amrex::Real b) noexcept
{
    using Real = amrex::Real;
    // Abscissae (15 points)
    static constexpr Real x[15] = {
        0.0000000000000000_rt,
       -0.9914553711208126_rt,
       -0.9491079123427585_rt,
       -0.8648644233597691_rt,
       -0.7415311855993945_rt,
       -0.5860872354676911_rt,
       -0.4058451513773972_rt,
       -0.2077849550078985_rt,
        0.2077849550078985_rt,
        0.4058451513773972_rt,
        0.5860872354676911_rt,
        0.7415311855993945_rt,
        0.8648644233597691_rt,
        0.9491079123427585_rt,
        0.9914553711208126_rt
    };
    // Kronrod weights
    static constexpr Real wK[15] = {
        0.2094821410847278_rt,
        0.02293532201052922_rt,
        0.06309209262997855_rt,
        0.1047900103222502_rt,
        0.1406532597155259_rt,
        0.1690047266392679_rt,
        0.1903505780647854_rt,
        0.2044329400752989_rt,
        0.2044329400752989_rt,
        0.1903505780647854_rt,
        0.1690047266392679_rt,
        0.1406532597155259_rt,
        0.1047900103222502_rt,
        0.06309209262997855_rt,
        0.02293532201052922_rt
    };
    // Gauss weights embedded in Kronrod nodes
    static constexpr Real wG[15] = {
        0.4179591836734694_rt, // x=0
        0.0000000000000000_rt,
        0.1294849661688697_rt,
        0.0000000000000000_rt,
        0.2797053914892766_rt,
        0.0000000000000000_rt,
        0.3818300505051189_rt,
        0.0000000000000000_rt,
        0.0000000000000000_rt,
        0.3818300505051189_rt,
        0.0000000000000000_rt,
        0.2797053914892766_rt,
        0.0000000000000000_rt,
        0.1294849661688697_rt,
        0.0000000000000000_rt
    };

    Real c = (a + b) * 0.5_rt;
    Real h = (b - a) * 0.5_rt;
    Real sumK = 0.0_rt;
    Real sumG = 0.0_rt;
    for (int i = 0; i < 15; ++i) {
        Real xi = x[i];
        Real f1 = f(c + h * xi);
        Real fsum = f1;
        if (xi != static_cast<Real>(0.0_rt)) {
            fsum += f(c - h * xi);
        }
        sumK += wK[i] * fsum;
        sumG += wG[i] * fsum;
    }
    Real integralK = sumK * h;
    Real integralG = sumG * h;
    Real error = amrex::Math::abs(integralK - integralG);
    return {integralK, error};
}

#endif // MICROPHYSICS_QUADRATURE_H