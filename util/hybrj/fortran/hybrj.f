      subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode,
     *                 factor,nprint,info,nfev,njev,r,lr,qtf,wa1,wa2,
     *                 wa3,wa4)
      integer n,ldfjac,maxfev,mode,nprint,info,nfev,njev,lr
      double precision xtol,factor
      double precision x(n),fvec(n),fjac(ldfjac,n),diag(n),r(lr),
     *                 qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
c     **********
c
c     subroutine hybrj
c
c     the purpose of hybrj is to find a zero of a system of
c     n nonlinear functions in n variables by a modification
c     of the powell hybrid method. the user must provide a
c     subroutine which calculates the functions and the jacobian.
c
c     the subroutine statement is
c
c       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,
c                        mode,factor,nprint,info,nfev,njev,r,lr,qtf,
c                        wa1,wa2,wa3,wa4)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions and the jacobian. fcn must
c         be declared in an external statement in the user
c         calling program, and should be written as follows.
c
c         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag)
c         integer n,ldfjac,iflag
c         double precision x(n),fvec(n),fjac(ldfjac,n)
c         ----------
c         if iflag = 1 calculate the functions at x and
c         return this vector in fvec. do not alter fjac.
c         if iflag = 2 calculate the jacobian at x and
c         return this matrix in fjac. do not alter fvec.
c         ---------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of hybrj.
c         in this case set iflag to a negative integer.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... dogleg,dpmpar,enorm,
c                            qform,qrfac,r1mpyq,r1updt
c
c       fortran-supplied ... dabs,dmax1,dmin1,mod
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i,iflag,iter,j,jm1,l,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm,
     *                 prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,
     *                 zero
      double precision dpmpar,enorm
      data one,p1,p5,p001,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/

      logical test

      logical finished

      finished = .false.

c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      info = 0
      iflag = 0
      nfev = 0
      njev = 0
c
c     check the input parameters for errors.
c

      test = (n .le. 0 .or. ldfjac .lt. n .or. xtol .lt. zero
     *    .or. maxfev .le. 0 .or. factor .le. zero
     *    .or. lr .lt. (n*(n + 1))/2)

      if (.not. test) then

         if (mode == 2) then
            do j = 1, n
               if (diag(j) .le. zero) then
                  finished = .true.
                  exit
               end if
            end do
         end if
c
c     evaluate the function at the starting point
c     and calculate its norm.
c
         iflag = 1
         call fcn(n,x,fvec,fjac,ldfjac,iflag)
         nfev = 1
         if (iflag .lt. 0) then
            finished = .true.
         else
            fnorm = enorm(n,fvec)
         end if
c
c     initialize iteration counter and monitors.
c
         iter = 1
         ncsuc = 0
         ncfail = 0
         nslow1 = 0
         nslow2 = 0
c
c     beginning of the outer loop.
c
         do while (.not. finished)
            jeval = .true.
c
c        calculate the jacobian matrix.
c
            iflag = 2
            call fcn(n,x,fvec,fjac,ldfjac,iflag)
            njev = njev + 1
            if (iflag .lt. 0) then
               finished = .true.
               exit
            endif
c
c        compute the qr factorization of the jacobian.
c
            call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
c
c        on the first iteration and if mode is 1, scale according
c        to the norms of the columns of the initial jacobian.
c
            if (iter == 1) then
               if (mode /= 2) then
                  do j = 1, n
                     diag(j) = wa2(j)
                     if (wa2(j) .eq. zero) diag(j) = one
                  end do
               end if
c
c        on the first iteration, calculate the norm of the scaled x
c        and initialize the step bound delta.
c
               do j = 1, n
                  wa3(j) = diag(j)*x(j)
               end do
               xnorm = enorm(n,wa3)
               delta = factor*xnorm
               if (delta .eq. zero) delta = factor
            end if
c
c        form (q transpose)*fvec and store in qtf.
c
            do i = 1, n
               qtf(i) = fvec(i)
            end do
            do j = 1, n
               if (fjac(j,j) /= zero) then
                  sum = zero
                  do i = j, n
                     sum = sum + fjac(i,j)*qtf(i)
                  end do
                  temp = -sum/fjac(j,j)
                  do i = j, n
                     qtf(i) = qtf(i) + fjac(i,j)*temp
                  end do
               end if
            end do
c
c        copy the triangular factor of the qr factorization into r.
c
            sing = .false.
            do j = 1, n
               l = j
               jm1 = j - 1
               if (jm1 >= 1) then
                  do i = 1, jm1
                     r(l) = fjac(i,j)
                     l = l + n - i
                  end do
               end if
               r(l) = wa1(j)
               if (wa1(j) .eq. zero) sing = .true.
            end do
c
c        accumulate the orthogonal factor in fjac.
c
            call qform(n,n,fjac,ldfjac,wa1)
c
c        rescale if necessary.
c
            if (mode /= 2) then
               do j = 1, n
                  diag(j) = dmax1(diag(j),wa2(j))
               end do
            end if
c
c        beginning of the inner loop.
c
            do while (.true.)
c
c           if requested, call fcn to enable printing of iterates.
c
               if (nprint > 0) then
                  iflag = 0
                  if (mod(iter-1,nprint) .eq. 0)
     *                 call fcn(n,x,fvec,fjac,ldfjac,iflag)
                  if (iflag .lt. 0) then
                     finished = .true.
                     exit
                  endif
               end if
c
c           determine the direction p.
c
               call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
c
c           store the direction p and x + p. calculate the norm of p.
c
               do j = 1, n
                  wa1(j) = -wa1(j)
                  wa2(j) = x(j) + wa1(j)
                  wa3(j) = diag(j)*wa1(j)
               end do
               pnorm = enorm(n,wa3)
c
c           on the first iteration, adjust the initial step bound.
c
               if (iter .eq. 1) delta = dmin1(delta,pnorm)
c
c           evaluate the function at x + p and calculate its norm.
c
               iflag = 1
               call fcn(n,wa2,wa4,fjac,ldfjac,iflag)
               nfev = nfev + 1
               if (iflag .lt. 0) then
                  finished = .true.
                  exit
               end if
               fnorm1 = enorm(n,wa4)
c     
c           compute the scaled actual reduction.
c
               actred = -one
               if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
c           compute the scaled predicted reduction.
c
               l = 1
               do i = 1, n
                  sum = zero
                  do j = i, n
                     sum = sum + r(l)*wa1(j)
                     l = l + 1
                  end do
                  wa3(i) = qtf(i) + sum
               end do
               temp = enorm(n,wa3)
               prered = zero
               if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
c
c           compute the ratio of the actual to the predicted
c           reduction.
c
               ratio = zero
               if (prered .gt. zero) ratio = actred/prered
c
c           update the step bound.
c
               if (ratio < p1) then
                  ncsuc = 0
                  ncfail = ncfail + 1
                  delta = p5*delta
               else
                  ncfail = 0
                  ncsuc = ncsuc + 1
                  if (ratio .ge. p5 .or. ncsuc .gt. 1)
     *                 delta = dmax1(delta,pnorm/p5)
                  if (dabs(ratio-one) .le. p1) delta = pnorm/p5
               end if
c
c           test for successful iteration.
c
               if (ratio >= p0001) then
c
c           successful iteration. update x, fvec, and their norms.
c
                  do j = 1, n
                     x(j) = wa2(j)
                     wa2(j) = diag(j)*x(j)
                     fvec(j) = wa4(j)
                  end do
                  xnorm = enorm(n,wa2)
                  fnorm = fnorm1
                  iter = iter + 1
               end if
c
c           determine the progress of the iteration.
c
               nslow1 = nslow1 + 1
               if (actred .ge. p001) nslow1 = 0
               if (jeval) nslow2 = nslow2 + 1
               if (actred .ge. p1) nslow2 = 0
c     
c           test for convergence.
c
               if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
               if (info .ne. 0) then
                  finished = .true.
                  exit
               end if
c
c           tests for termination and stringent tolerances.
c
               if (nfev .ge. maxfev) info = 2
               if (p1*dmax1(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
               if (nslow2 .eq. 5) info = 4
               if (nslow1 .eq. 10) info = 5
               if (info .ne. 0) then
                  finished = .true.
                  exit
               end if
c
c           criterion for recalculating jacobian.
c
               if (ncfail .eq. 2) exit
c
c           calculate the rank one modification to the jacobian
c           and update qtf if necessary.
c
               do j = 1, n
                  sum = zero
                  do i = 1, n
                     sum = sum + fjac(i,j)*wa4(i)
                  end do
                  wa2(j) = (sum - wa3(j))/pnorm
                  wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
                  if (ratio .ge. p0001) qtf(j) = sum
               end do
c
c           compute the qr factorization of the updated jacobian.
c
               call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
               call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
               call r1mpyq(1,n,qtf,1,wa2,wa3)
c
c           end of the inner loop.
c
               jeval = .false.
            end do
         
c
c        end of the outer loop.
c

            if (finished) then
               exit
            end if
         end do
      endif

c
c     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(n,x,fvec,fjac,ldfjac,iflag)
      return
c
c     last card of subroutine hybrj.
c
      end
