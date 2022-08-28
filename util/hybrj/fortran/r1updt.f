      subroutine r1updt(m,n,s,ls,u,v,w,sing)
      integer m,n,ls
      logical sing
      double precision s(ls),u(m),v(n),w(m)
c     **********
c
c     subroutine r1updt
c
c     given an m by n lower trapezoidal matrix s, an m-vector u,
c     and an n-vector v, the problem is to determine an
c     orthogonal matrix q such that
c
c                   t
c           (s + u*v )*q
c
c     is again lower trapezoidal.
c
c     this subroutine determines q as the product of 2*(n - 1)
c     transformations
c
c           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
c     where gv(i), gw(i) are givens rotations in the (i,n) plane
c     which eliminate elements in the i-th and n-th planes,
c     respectively. q itself is not accumulated, rather the
c     information to recover the gv, gw rotations is returned.
c
c     the subroutine statement is
c
c       subroutine r1updt(m,n,s,ls,u,v,w,sing)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of s.
c
c       n is a positive integer input variable set to the number
c         of columns of s. n must not exceed m.
c
c       s is an array of length ls. on input s must contain the lower
c         trapezoidal matrix s stored by columns. on output s contains
c         the lower trapezoidal matrix produced as described above.
c
c       ls is a positive integer input variable not less than
c         (n*(2*m-n+1))/2.
c
c       u is an input array of length m which must contain the
c         vector u.
c
c       v is an array of length n. on input v must contain the vector
c         v. on output v(i) contains the information necessary to
c         recover the givens rotation gv(i) described above.
c
c       w is an output array of length m. w(i) contains information
c         necessary to recover the givens rotation gw(i) described
c         above.
c
c       sing is a logical output variable. sing is set true if any
c         of the diagonal elements of the output s are zero. otherwise
c         sing is set false.
c
c     subprograms called
c
c       minpack-supplied ... dpmpar
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more,
c     john l. nazareth
c
c     **********
      integer i,j,jj,l,nmj,nm1
      double precision cos,cotan,giant,one,p5,p25,sin,tan,tau,temp,
     *                 zero
      double precision dpmpar
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
c
c     giant is the largest magnitude.
c
      giant = dpmpar(3)
c
c     initialize the diagonal element pointer.
c
      jj = (n*(2*m - n + 1))/2 - (m - n)
c
c     move the nontrivial part of the last column of s into w.
c
      l = jj
      do i = n, m
         w(i) = s(l)
         l = l + 1
      end do
c
c     rotate the vector v into a multiple of the n-th unit vector
c     in such a way that a spike is introduced into w.
c
      nm1 = n - 1
      if (nm1 >= 1) then
         do nmj = 1, nm1
            j = n - nmj
            jj = jj - (m - j + 1)
            w(j) = zero
            if (v(j) /= zero) then
c
c        determine a givens rotation which eliminates the
c        j-th element of v.
c
               if (dabs(v(n)) < dabs(v(j))) then
                  cotan = v(n)/v(j)
                  sin = p5/dsqrt(p25+p25*cotan**2)
                  cos = sin*cotan
                  tau = one
                  if (dabs(cos)*giant .gt. one) tau = one/cos
               else
                  tan = v(j)/v(n)
                  cos = p5/dsqrt(p25+p25*tan**2)
                  sin = cos*tan
                  tau = sin
               end if
c
c        apply the transformation to v and store the information
c        necessary to recover the givens rotation.
c
               v(n) = sin*v(j) + cos*v(n)
               v(j) = tau
c
c        apply the transformation to s and extend the spike in w.
c
               l = jj
               do i = j, m
                  temp = cos*s(l) - sin*w(i)
                  w(i) = sin*s(l) + cos*w(i)
                  s(l) = temp
                  l = l + 1
               end do
            end if
         end do
      end if
c
c     add the spike from the rank 1 update to w.
c
      do i = 1, m
         w(i) = w(i) + v(n)*u(i)
      end do
c
c     eliminate the spike.
c
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
