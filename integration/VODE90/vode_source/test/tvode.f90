program tvode

  use amrex_fort_module, only : rt => amrex_real
  use dvode_module, only: dvode
  use dvode_type_module, only: dvode_t
  use tvode_rhs_module, only: FEX, JEX
  
  implicit none

  !  The following output was obtained from this program on a
  !  Cray-1 computer with the CFT compiler.
  ! 
  !  At t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
  !  At t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
  !  At t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
  !  At t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
  !  At t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
  !  At t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
  !  At t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
  !  At t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
  !  At t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
  !  At t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
  !  At t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
  !  At t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
  ! 
  !  No. steps = 595   No. f-s = 832   No. J-s =  13   No. LU-s = 112
  !   No. nonlinear iterations = 831
  !   No. nonlinear convergence failures =   0
  !   No. error test failures =  22
  
  type(dvode_t) :: dvode_state
  real(rt) :: RPAR(1), RTOL(1), T, TOUT
  real(rt) :: Y(3), ATOL(3), RWORK(67)
  integer    :: IWORK(33), NEQ, ITOL, ITASK, ISTATE, IOPT
  integer    :: LRW, LIW, MF, IOUT, IPAR(1)
  RPAR = 0.0d0
  IPAR = 0
  RWORK = 0.0d0
  IWORK = 0
  NEQ = 3
  Y(1) = 1.0D0
  Y(2) = 0.0D0
  Y(3) = 0.0D0
  T = 0.0D0
  TOUT = 0.4D0
  ITOL = 2
  RTOL(1) = 1.D-4
  ATOL(1) = 1.D-8
  ATOL(2) = 1.D-14
  ATOL(3) = 1.D-6
  ITASK = 1
  ISTATE = 1
  IOPT = 0
  LRW = 67
  LIW = 33
  MF = 21
  do IOUT = 1,12
     ISTATE = 1
     CALL DVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE, &
          IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR, dvode_state)
     WRITE(6,20)T,Y(1),Y(2),Y(3)
20   FORMAT(' At t =',D12.4,'   y =',3D14.6)
     IF (ISTATE .LT. 0) GO TO 80
     TOUT = TOUT*10.
  end do
  WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19), &
       IWORK(20),IWORK(21),IWORK(22)
60 FORMAT(/' No. steps =',I4,'   No. f-s =',I4, &
        '   No. J-s =',I4,'   No. LU-s =',I4/ &
        '  No. nonlinear iterations =',I4/ &
        '  No. nonlinear convergence failures =',I4/ &
        '  No. error test failures =',I4/)
  STOP
80 WRITE(6,90)ISTATE
90 FORMAT(///' Error halt: ISTATE =',I3)
  STOP
end program tvode
 
