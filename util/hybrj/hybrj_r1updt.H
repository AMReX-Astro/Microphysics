#ifndef HYBRJ_R1UPDT_H
#define HYBRJ_R1UPDT_H

#include <hybrj_type.H>

// given an n by n lower trapezoidal matrix s, an n-vector u,
// and an n-vector v, the problem is to determine an
// orthogonal matrix q such that
//
//                   t
//           (s + u*v )*q
//
// is again lower trapezoidal.
//
// this subroutine determines q as the product of 2*(n - 1)
// transformations
//
//          gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
//
// where gv(i), gw(i) are givens rotations in the (i,n) plane
// which eliminate elements in the i-th and n-th planes,
// respectively. q itself is not accumulated, rather the
// information to recover the gv, gw rotations is returned.

// changes:
// we now assume that the matrix is nxn
//
// s is an array of length (n*(n+1))/2. on input s must contain the lower
//   trapezoidal matrix s stored by columns. on output s contains
//   the lower trapezoidal matrix produced as described above.
//
// u is an input array of length m which must contain the
//   vector u.
//
// v is an array of length n. on input v must contain the vector
//   v. on output v(i) contains the information necessary to
//   recover the givens rotation gv(i) described above.
//
//  w is an output array of length m. w(i) contains information
//    necessary to recover the givens rotation gw(i) described
//    above.
//
//  sing is a logical output variable. sing is set true if any
//       of the diagonal elements of the output s are zero. otherwise
//       sing is set false.
//
// argonne national laboratory. minpack project. march 1980.
// burton s. garbow, kenneth e. hillstrom, jorge j. more,
// john l. nazareth

template<int neqs>
void r1updt(amrex::Array1D<Real, 1, (neqs * (neqs + 1)) / 2>& s,
            amrex::Array1D<Real, 1, neqs>& u,
            amrex::Array1D<Real, 1, neqs>& v,
            amrex::Array1D<Real, 1, neqs>& w,
            bool& sing) {

    const Real giant = std::numeric_limits<Real>::max();

    // initialize the diagonal element pointer.

    int jj = (neqs * (neqs + 1))/2;

    // move the nontrivial part of the last column of s into w.

    int l = jj;
    for (int i = neqs; i <= neqs; ++i) {
         w(i) = s(l);
         l++;
    }

    // rotate the vector v into a multiple of the n-th unit vector
    // in such a way that a spike is introduced into w.

    int nm1 = neqs - 1;
    if (nm1 >= 1) {
        for (int nmj = 1; nmj <= nm1; ++nmj) {
            int j = neqs - nmj;
            jj += -(neqs - j + 1);
            w(j) = 0.0_rt;
            if (v(j) != 0.0_rt) {

                // determine a givens rotation which eliminates the
                // j-th element of v.

                Real tau;
                Real fcos;
                Real fsin;

                if (std::abs(v(neqs)) < std::abs(v(j))) {
                    Real cotan = v(neqs) / v(j);
                    fsin = 0.5_rt / std::sqrt(0.25_rt + 0.25_rt * cotan * cotan);
                    fcos = fsin * cotan;
                    tau = 1.0_rt;
                    if (std::abs(fcos) * giant > 1.0_rt) {
                        tau = 1.0_rt / fcos;
                    }
                } else {
                    Real ftan = v(j) / v(neqs);
                    fcos = 0.5_rt / std::sqrt(0.25_rt + 0.25_rt * ftan * ftan);
                    fsin = fcos * ftan;
                    tau = fsin;
                }

                // apply the transformation to v and store the information
                // necessary to recover the givens rotation.

                v(neqs) = fsin * v(j) + fcos * v(neqs);
                v(j) = tau;

                // apply the transformation to s and extend the spike in w.

                l = jj;
                for (int i = j; i <= neqs; ++i) {
                    Real temp = fcos * s(l) - fsin * w(i);
                    w(i) = fsin * s(l) + fcos * w(i);
                    s(l) = temp;
                    l++;
                }
            }
        }
    }

    // add the spike from the rank 1 update to w.

    for (int i = 1; i <= neqs; ++i) {
        w(i) += v(neqs) * u(i);
    }

    // eliminate the spike.

    sing = false;

    if (nm1 >= 1) {
        for (int j = 1; j <= nm1; ++j) {
            if (w(j) != 0.0_rt) {

                // determine a givens rotation which eliminates the
                // j-th element of the spike.

                Real tau;
                Real fcos;
                Real fsin;

                if (std::abs(s(jj)) < std::abs(w(j))) {
                    Real cotan = s(jj) / w(j);
                    fsin = 0.5_rt / std::sqrt(0.25_rt + 0.25_rt * cotan * cotan);
                    fcos = fsin * cotan;
                    tau = 1.0_rt;
                    if (std::abs(fcos) * giant > 1.0_rt) {
                        tau = 1.0_rt / fcos;
                    }
                } else {
                    Real ftan = w(j) / s(jj);
                    fcos = 0.5_rt / std::sqrt(0.25_rt + 0.25_rt * ftan * ftan);
                    fsin = fcos * ftan;
                    tau = fsin;
                }

                // apply the transformation to s and reduce the spike in w.

                l = jj;
                for (int i = j; i <= neqs; ++i) {
                    Real temp = fcos * s(l) + fsin * w(i);
                    w(i) = -fsin * s(l) + fcos * w(i);
                    s(l) = temp;
                    l++;
                }

                // store the information necessary to recover the
                // givens rotation.

                w(j) = tau;
            }

            // test for zero diagonal elements in the output s.

            if (s(jj) == 0.0_rt) {
                sing = true;
            }
            jj += (neqs - j + 1);
        }
    }

    // move w back into the last column of the output s.

    l = jj;
    for (int i = neqs; i <= neqs; ++i) {
        s(l) = w(i);
        l++;
    }

    if (s(jj) == 0.0_rt) {
        sing = true;
    }
}
#endif
