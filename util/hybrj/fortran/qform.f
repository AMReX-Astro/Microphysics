      subroutine qform(m,n,q,ldq,wa)
      integer m,n,ldq
      double precision q(ldq,m),wa(m)
c     **********
c
c     subroutine qform
c
c     this subroutine proceeds from the computed qr factorization of
c     an m by n matrix a to accumulate the m by m orthogonal matrix
c     q from its factored form.
c
c     the subroutine statement is
c
c       subroutine qform(m,n,q,ldq,wa)
c
c     where
c
c       m is a positive integer input variable set to the number
c         of rows of a and the order of q.
c
c       n is a positive integer input variable set to the number
c         of columns of a.
c
c       q is an m by m array. on input the full lower trapezoid in
c         the first min(m,n) columns of q contains the factored form.
c         on output q has been accumulated into a square matrix.
c
c       ldq is a positive integer input variable not less than m
c         which specifies the leading dimension of the array q.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       fortran-supplied ... min0
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,j,jm1,k,l,minmn,np1
      double precision one,sum,temp,zero
      data one,zero /1.0d0,0.0d0/
c
c     zero out upper triangle of q in the first min(m,n) columns.
c
      minmn = min0(m,n)
      if (minmn >= 2) then
         do j = 2, minmn
            jm1 = j - 1
            do i = 1, jm1
               q(i,j) = zero
            end do
         end do
      end if
c
c     initialize remaining columns to those of the identity matrix.
c
      np1 = n + 1
      if (m >= np1) then
         do j = np1, m
            do i = 1, m
               q(i,j) = zero
            end do
            q(j,j) = one
         end do
      end if
c
c     accumulate q from its factored form.
c
      do l = 1, minmn
         k = minmn - l + 1
         do i = k, m
            wa(i) = q(i,k)
            q(i,k) = zero
         end do
         q(k,k) = one
         if (wa(k) /= zero) then
            do j = k, m
               sum = zero
               do i = k, m
                  sum = sum + q(i,j)*wa(i)
               end do
               temp = sum/wa(k)
               do i = k, m
                  q(i,j) = q(i,j) - temp*wa(i)
               end do
            end do
         end if
      end do
      return
c
c     last card of subroutine qform.
c
      end
