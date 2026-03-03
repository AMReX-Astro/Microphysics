program dvode_driver
  implicit none
  external :: fex, jex
  integer, parameter :: neq = 3
  integer :: itol, itask, istate, iopt, lrw, liw, mf, iout
  integer :: iwork(1000), ipar(1)
  double precision :: y(neq), t, tout, rtol, atol(neq), rwork(2000), rpar(1)

  y(1) = 1.0d0
  y(2) = 0.0d0
  y(3) = 0.0d0
  t = 0.0d0
  tout = 0.4d0

  itol = 2
  rtol = 1.0d-8
  atol(1) = 1.0d-14
  atol(2) = 1.0d-14
  atol(3) = 1.0d-14

  itask = 1
  istate = 1
  iopt = 0
  lrw = 2000
  liw = 1000
  mf = 21

  rwork = 0.0d0
  iwork = 0

  do iout = 1, 12
     call dvode(fex, neq, y, t, tout, itol, rtol, atol, itask, istate, &
                iopt, rwork, lrw, iwork, liw, jex, mf, rpar, ipar)
     write(6, '(A,1PE20.12,3(1PE20.12))') '[DRIVER] t,y=', t, y(1), y(2), y(3)
     if (istate < 0) then
        write(6, '(A,I6)') '[DRIVER] ISTATE=', istate
        stop 1
     end if
     tout = tout * 10.0d0
  end do

  write(6, '(A,3(1PE20.12))') '[DRIVER_FINAL] y=', y(1), y(2), y(3)
  write(6, '(A,3(I10))') '[DRIVER_STATS] NST NFE NJE=', iwork(11), iwork(12), iwork(13)
end program dvode_driver

subroutine fex (neq, t, y, ydot, rpar, ipar)
  implicit none
  integer neq, ipar(*)
  double precision t, y(neq), ydot(neq), rpar(*)
  ydot(1) = -0.04d0 * y(1) + 1.0d4 * y(2) * y(3)
  ydot(3) = 3.0d7 * y(2) * y(2)
  ydot(2) = -ydot(1) - ydot(3)
end subroutine fex

subroutine jex (neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
  implicit none
  integer neq, ml, mu, nrpd, ipar(*)
  double precision t, y(neq), pd(nrpd,neq), rpar(*)

  pd(1,1) = -0.04d0
  pd(1,2) = 1.0d4 * y(3)
  pd(1,3) = 1.0d4 * y(2)

  pd(2,1) = 0.04d0
  pd(2,2) = -1.0d4 * y(3) - 6.0d7 * y(2)
  pd(2,3) = -1.0d4 * y(2)

  pd(3,1) = 0.0d0
  pd(3,2) = 6.0d7 * y(2)
  pd(3,3) = 0.0d0
end subroutine jex
