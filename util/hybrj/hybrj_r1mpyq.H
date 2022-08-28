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
// * redone to require an nxn matrix

template<int neqs>
void (amrex::Array2D<Real, 1, neqs, 1, neqs>& a,
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

    Real cos;
    Real sin;

    if (nm1 >= 1) {
        for (int nmj = 1; nmj <= nm1; ++nmj) {
            int j = n - nmj;
            if (std::abs(v(j)) > 1.0_rt) {
                cos = 1.0_rt / v(j);
                sin = std::sqrt(1.0_rt - cos*cos);
            }
            if (std::abs(v(j)) <= 1.0_rt) {
                sin = v(j);
                cos = std::sqrt(1.0_rt - sin*sin);
            }
            for (int i = 1; i <= neqs; ++i) {
                Real temp = cos * a(i,j) - sin * a(i,n);
                a(i,n) = sin * a(i,j) + cos * a(i,n);
                a(i,j) = temp;
            }
        }

        // apply the second set of givens rotations to a.

        for (int j = 1; j <= nm1; ++j) {
            if (std::abs(w(j)) > 1.0_rt) {
                cos = 1.0_rt / w(j);
                sin = std::sqrt(1.0_rt - cos*cos);
            }
            if (std::abs(w(j)) <= 1.0_rt) {
                sin = w(j);
                cos = std::sqrt(1.0_rt - sin*sin);
            }
            for (int i = 1; i <= neqs; ++i) {
                Real temp = cos * a(i,j) + sin * a(i,n);
                a(i,n) = -sin*a(i,j) + cos*a(i,n);
                a(i,j) = temp;
            }
        }
    }

}
#endif
