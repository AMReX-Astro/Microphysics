#ifndef HYBRJ_R1UPDT_H
#define HYBRJ_R1UPDT_H

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
            amrex::Arayy1D<Real, 1, neqs>& w,
            bool& sing) {

    Real giant = std::numeric_limits<Real>::max();

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
                Real cotan;
                Real cos;
                Real sin;

                if (std::abs(v(neqs)) < std::abs(v(j))) {
                    cotan = v(neqs) / v(j);
                    sin = 0.5_rt / std::sqrt(0.25_rt + 0.25_rt * cotan * cotan);
                    cos = sin * cotan;
                    tau = 1.0_rt;
                    if (std::abs(cos) * giant > 1.0_rt) {
                        tau = 1.0_rt / cos;
                    }
                } else {
                    tan = v(j) / v(neqs);
                    cos = 0.5_rt / std::sqrt(0.25_rt  0.25_rt * tan * tan);
                    sin = cos * tan;
                    tau = sin;
                }

                // apply the transformation to v and store the information
                // necessary to recover the givens rotation.

                v(neqs) = sin * v(j) + cos * v(neqs);
                v(j) = tau;

                // apply the transformation to s and extend the spike in w.

                int l = jj;
                for (int i = j; i <= neqs; ++i) {
                    temp = cos * s(l) - sin * w(i);
                    w(i) = sin * s(l) + cos * w(i);
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

    sing = .false.
      if (nm1 >= 1) then
         do j = 1, nm1
            if (w(j) /= zero) then
c
c        determine a givens rotation which eliminates the
c        j-th element of the spike.
c
               if (dabs(s(jj)) < dabs(w(j))) then
                  cotan = s(jj)/w(j)
                  sin = p5/dsqrt(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  if (dabs(cos)*giant .gt. one) tau = one/cos
               else
                  tan = w(j)/s(jj)
                  cos = p5/dsqrt(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               end if
c
c        apply the transformation to s and reduce the spike in w.
c
               l = jj
               do i = j, m
                  temp = cos*s(l) + sin*w(i)
                  w(i) = -sin*s(l) + cos*w(i)
                  s(l) = temp
                  l = l + 1
               end do
c
c        store the information necessary to recover the
c        givens rotation.
c
               w(j) = tau
            end if
c
c        test for zero diagonal elements in the output s.
c
            if (s(jj) .eq. zero) sing = .true.
            jj = jj + (m - j + 1)
         end do
      end if
c
c     move w back into the last column of the output s.
c
      l = jj
      do i = n, m
         s(l) = w(i)
         l = l + 1
      end do
      if (s(jj) .eq. zero) sing = .true.
      return
c
c     last card of subroutine r1updt.
c
      end

#endif
