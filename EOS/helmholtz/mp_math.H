#ifndef MP_MATH_H
#define MP_MATH_H

#include <AMReX_REAL.H>

namespace mp {

    // these are approximate functions based on the singularity EOS /
    // Not-Quite-Transcendental Functions and their Applications
    // by Miller et al.

    // approximation to log_2(x)
    inline auto fastlg2(const amrex::Real x) -> amrex::Real {
        // split into mant + exp
        int n{};
        const amrex::Real m = std::frexp(x, &n);
        return 2.0_rt * (m - 1.0_rt) + static_cast<amrex::Real>(n);
    }

    // approximation to 2**x
    inline auto fastpow2(const amrex::Real x) -> amrex::Real {
        // create a mantissa and exp compatible with our log
        const int flr = std::floor(x);
        const amrex::Real rem = x - flr;
        const amrex::Real mant = 0.5_rt * (rem + 1.0_rt);
        const int exponent = flr + 1;
        return std::ldexp(mant, exponent);
    }

}

#endif
