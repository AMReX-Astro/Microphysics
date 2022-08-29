!       The problem is to determine the values of x(1), x(2), ..., x(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*x(1))*x(1)           -2*x(2)                   = -1
!               -x(i-1) + (3-2*x(i))*x(i)         -2*x(i+1) = -1, i=2-8
!                                   -x(8) + (3-2*x(9))*x(9) = -1

!     DRIVER FOR HYBRJ EXAMPLE.

program test_hybrj

  INTEGER J,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,LR,NWRITE
  DOUBLE PRECISION XTOL,FACTOR,FNORM
  DOUBLE PRECISION X(9),FVEC(9),FJAC(9,9),DIAG(9),R(45),QTF(9)
  double precision WA1(9),WA2(9),WA3(9),WA4(9)
  DOUBLE PRECISION ENORM,DPMPAR
  EXTERNAL FCN

  N = 9

 !     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.

  DO J = 1, 9
     X(J) = -1.D0
  end do

  LDFJAC = 9
  LR = 45

  ! SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
  ! UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
  ! THIS IS THE RECOMMENDED SETTING.

  XTOL = DSQRT(DPMPAR(1))

  MAXFEV = 1000
  MODE = 2
  DO J = 1, 9
     DIAG(J) = 1.D0
  end do
  FACTOR = 1.D2
  NPRINT = 0

  CALL HYBRJ(FCN,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,DIAG, &
             MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,R,LR,QTF, &
             WA1,WA2,WA3,WA4)

  FNORM = ENORM(N,FVEC)
  print *, "final L2 norm of the residuals", FNORM
  print *, "number of function evaluations", NFEV
  print *, "number of Jacobian evaluations", NJEV
  print *, "exit parameter", INFO
  print *, "final approximate solution", (X(J),J=1,N)

END program test_hybrj


SUBROUTINE FCN(N,X,FVEC,FJAC,LDFJAC,IFLAG)
  INTEGER N,LDFJAC,IFLAG
  DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)

  ! SUBROUTINE FCN FOR HYBRJ EXAMPLE.

  INTEGER J,K
  DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
  DATA ZERO,ONE,TWO,THREE,FOUR /0.D0,1.D0,2.D0,3.D0,4.D0/

  IF (IFLAG == 0) return

  ! INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.

  IF (IFLAG /= 2) then
     DO K = 1, N
        TEMP = (THREE - TWO*X(K))*X(K)
        TEMP1 = ZERO
        IF (K .NE. 1) TEMP1 = X(K-1)
        TEMP2 = ZERO
        IF (K .NE. N) TEMP2 = X(K+1)
        FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
     enddo
  else
     DO K = 1, N
        DO J = 1, N
           FJAC(K,J) = ZERO
        end do
        FJAC(K,K) = THREE - FOUR*X(K)
        IF (K .NE. 1) FJAC(K,K-1) = -ONE
        IF (K .NE. N) FJAC(K,K+1) = -TWO
     end do
  end if

  RETURN
END SUBROUTINE FCN
