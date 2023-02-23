#ifndef HYBRJ_ENORM_H
#define HYBRJ_ENORM_H

#include <hybrj_type.H>
#include <AMReX_GpuUtility.H>

// given an n-vector x, this function calculates the
// euclidean norm of x.
//
// the euclidean norm is computed by accumulating the sum of
// squares in three different sums. the sums of squares for the
// small and large components are scaled so that no overflows
// occur. non-destructive underflows are permitted. underflows
// and overflows do not occur in the computation of the unscaled
// sum of squares for the intermediate components.
// the definitions of small, intermediate and large components
// depend on two constants, rdwarf and rgiant. the main
// restrictions on these constants are that rdwarf**2 not
// underflow and rgiant**2 not overflow. the constants
// given here are suitable for every known computer.
//
// argonne national laboratory. minpack project. march 1980.
// burton s. garbow, kenneth e. hillstrom, jorge j. more

template<int neqs>
Real enorm(const int n, amrex::Array1D<Real, 1, neqs>& x) {

    // note: n <= neqs

    const Real rdwarf = 3.834e-20_rt;
    const Real rgiant = 1.304e19_rt;

    Real s1 = 0.0_rt;
    Real s2 = 0.0_rt;
    Real s3 = 0.0_rt;

    Real x1max = 0.0_rt;
    Real x3max = 0.0_rt;

    Real agiant = rgiant / static_cast<Real>(n);

    for (int i = 1; i <= n; ++i) {
        Real xabs = std::abs(x(i));

	if ( amrex::isnan(xabs) || amrex::isinf(xabs) ) {
	  return 0.0_rt;
	}
        if (xabs <= rdwarf || xabs >= agiant) {
            if (xabs > rdwarf) {

                // sum for large components.

                if (xabs > x1max) {
                    s1 = 1.0_rt + s1 * std::pow(x1max/xabs, 2.0_rt);
                    x1max = xabs;
                } else {
                    s1 += std::pow(xabs/x1max, 2.0_rt);
                }
            } else {

                // sum for small components.

                if (xabs > x3max) {
                    s3 = 1.0_rt + s3 * std::pow(x3max/xabs, 2.0_rt);
                    x3max = xabs;
                } else {
                    if (xabs != 0.0_rt) {
                        s3 += std::pow(xabs/x3max, 2.0_rt);
                    }
                }
            }
        } else {

            // sum for intermediate components.

            s2 += xabs * xabs;
        }
    }

    // calculation of norm.

    Real _enorm = 0.0_rt;
    if (s1 != 0.0_rt) {
        _enorm = x1max * std::sqrt(s1 + (s2/x1max) / x1max);
    } else {
        if (s2 != 0.0_rt) {
            if (s2 >= x3max) {
                _enorm = std::sqrt(s2 * (1.0_rt + (x3max/s2)*(x3max*s3)));
            }
            if (s2 < x3max) {
                _enorm = std::sqrt(x3max * ((s2/x3max) + (x3max*s3)));
            }
        } else {
            _enorm = x3max * std::sqrt(s3);
        }
    }

    return _enorm;
}
#endif
