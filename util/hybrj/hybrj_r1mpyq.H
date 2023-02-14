#ifndef HYBRJ_R1MPYQ_H
#define HYBRJ_R1MPYQ_H

#include <hybrj.H>

// given an m by n matrix a, this subroutine computes a*q where
//  q is the product of 2*(n - 1) transformations
//
//       gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
//
// and gv(i), gw(i) are givens rotations in the (i,n) plane which
// eliminate elements in the i-th and n-th planes, respectively.
// q itself is not given, rather the information to recover the
// gv, gw rotations is supplied.
//
// changes:
// * Split into 2 versions, one takes an nxn matrix and one takes a vector of length n

template<int neqs>
void r1mpyq(amrex::Array2D<Real, 1, neqs, 1, neqs>& a,
            amrex::Array1D<Real, 1, neqs>& v,
            amrex::Array1D<Real, 1, neqs>& w) {

    // a is an neqs by neqs array. on input a must contain the matrix
    //   to be postmultiplied by the orthogonal matrix q
    //   described above. on output a*q has replaced a.
    //
    // v is an input array of length neqs. v(i) must contain the
    //   information necessary to recover the givens rotation gv(i)
    //   described above.
    //
    // w is an input array of length neqs. w(i) must contain the
    //   information necessary to recover the givens rotation gw(i)
    //   described above.
    //
    // argonne national laboratory. minpack project. march 1980.
    // burton s. garbow, kenneth e. hillstrom, jorge j. more

    // apply the first set of givens rotations to a.

    int nm1 = neqs - 1;

    Real fcos;
    Real fsin;

    if (nm1 >= 1) {
        for (int nmj = 1; nmj <= nm1; ++nmj) {
            int j = neqs - nmj;
            if (std::abs(v(j)) > 1.0_rt) {
                fcos = 1.0_rt / v(j);
                fsin = std::sqrt(1.0_rt - fcos*fcos);
            }
	    else {
            // if (std::abs(v(j)) <= 1.0_rt) {
                fsin = v(j);
                fcos = std::sqrt(1.0_rt - fsin*fsin);
            }
            for (int i = 1; i <= neqs; ++i) {
                Real temp = fcos * a(i,j) - fsin * a(i,neqs);
                a(i,neqs) = fsin * a(i,j) + fcos * a(i,neqs);
                a(i,j) = temp;
            }
        }

        // apply the second set of givens rotations to a.

        for (int j = 1; j <= nm1; ++j) {
            if (std::abs(w(j)) > 1.0_rt) {
                fcos = 1.0_rt / w(j);
                fsin = std::sqrt(1.0_rt - fcos*fcos);
            }
	    else {
            // if (std::abs(w(j)) <= 1.0_rt) {
                fsin = w(j);
                fcos = std::sqrt(1.0_rt - fsin*fsin);
            }
            for (int i = 1; i <= neqs; ++i) {
                Real temp = fcos * a(i,j) + fsin * a(i,neqs);
                a(i,neqs) = -fsin*a(i,j) + fcos*a(i,neqs);
                a(i,j) = temp;
            }
        }
    }

}


template<int neqs>
void r1mpyq(amrex::Array1D<Real, 1, neqs>& a,
            amrex::Array1D<Real, 1, neqs>& v,
            amrex::Array1D<Real, 1, neqs>& w) {

    // a is an neqs vector.  This is based on the original version
    // with m = 1.
    //
    // v is an input array of length neqs. v(i) must contain the
    //   information necessary to recover the givens rotation gv(i)
    //   described above.
    //
    // w is an input array of length neqs. w(i) must contain the
    //   information necessary to recover the givens rotation gw(i)
    //   described above.
    //
    // argonne national laboratory. minpack project. march 1980.
    // burton s. garbow, kenneth e. hillstrom, jorge j. more

    // apply the first set of givens rotations to a.


    int nm1 = neqs - 1;

    Real fcos;
    Real fsin;

    if (nm1 >= 1) {
        for (int nmj = 1; nmj <= nm1; ++nmj) {
            int j = neqs - nmj;
            if (std::abs(v(j)) > 1.0_rt) {
                fcos = 1.0_rt / v(j);
                fsin = std::sqrt(1.0_rt - fcos*fcos);
            }
	    else {
            // if (std::abs(v(j)) <= 1.0_rt) {
                fsin = v(j);
                fcos = std::sqrt(1.0_rt - fsin*fsin);
            }
            Real temp = fcos*a(j) - fsin*a(neqs);
            a(neqs) = fsin*a(j) + fcos*a(neqs);
            a(j) = temp;
        }

        // apply the second set of givens rotations to a.

        for (int j = 1; j <= nm1; ++j) {
            if (std::abs(w(j)) > 1.0_rt) {
                fcos = 1.0_rt / w(j);
                fsin = std::sqrt(1.0_rt - fcos*fcos);
            }
	    else {
            // if (std::abs(w(j)) <= 1.0_rt) {
                fsin = w(j);
                fcos = std::sqrt(1.0_rt - fsin*fsin);
            }
            Real temp = fcos*a(j) + fsin*a(neqs);
            a(neqs) = -fsin*a(j) + fcos*a(neqs);
            a(j) = temp;
        }
    }
}
#endif
