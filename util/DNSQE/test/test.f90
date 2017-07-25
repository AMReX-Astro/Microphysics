program test

  ! the test example from the dnsqe source comments

  integer :: j, n, iopt, nprint, info, lwa
  double precision :: tol, fnorm
  double precision :: x(9), fvec(9), wa(180)
  double precision :: denorm, d1mach

  external fcn

  iopt = 2
  n = 9

  ! the following starting values provide a rough solution.
  do j = 1, 9
     x(j) = -1.e0
  enddo

  lwa = 180
  nprint = 0

  ! set tol to the square root of the machine precision.  unless high
  ! precision solutions are required, this is the recommended setting.
  tol = sqrt(d1mach(4))

  call dnsqe(fcn, jac, iopt, n, x, fvec, tol, nprint, info, wa, lwa)
  fnorm = denorm(n,fvec)
  write (* ,1000) fnorm, info, (x(j),j=1,n)
  
1000 format (5x,' final l2 norm of the residuals',e15.7 // &
             5x,' exit parameter',16x,i10 // &
             5x,' final approximate solution' // (5x,3e15.7))
end program test

subroutine fcn(n, x, fvec, iflag)

  integer :: n, iflag
  double precision :: x(n), fvec(n)
  integer :: k
  double precision :: temp, temp1, temp2
  double precision, parameter :: zero = 0.0d0
  double precision, parameter :: one = 1.0d0
  double precision, parameter :: two = 2.0d0
  double precision, parameter :: three = 3.0d0

  do k = 1, n
     temp = (three - two*x(k))*x(k)
     temp1 = zero
     if (k .ne. 1) temp1 = x(k-1)
     temp2 = zero
     if (k .ne. n) temp2 = x(k+1)
     fvec(k) = temp - temp1 - two*temp2 + one
  enddo

end subroutine fcn
