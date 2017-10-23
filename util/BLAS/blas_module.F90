module blas_module

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif
  SUBROUTINE DCOPYN(N,DX,INCX,DY,INCY)

  ! Only operates on vectors of size N

    INTEGER INCX,INCY,N
    DOUBLE PRECISION DX(N),DY(N)
! *  Purpose
! *  =======
! *
! *     copies a vector, x, to a vector, y.
! *     uses unrolled loops for increments equal to one.
! *     jack dongarra, linpack, 3/11/78.
! *     modified 12/3/93, array(1) declarations changed to array(*)
    INTEGER I,IX,IY,M,MP1

    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20

    ! code for unequal increments or equal increments
    ! not equal to 1

    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
       DY(IY) = DX(IX)
       IX = IX + INCX
       IY = IY + INCY
    end do
    RETURN

    ! code for both increments equal to 1

    ! clean-up loop

20  M = MOD(N,7)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
       DY(I) = DX(I)
    end do
    IF (N.LT.7) RETURN
40  MP1 = M + 1
    DO I = MP1,N,7
       DY(I) = DX(I)
       DY(I+1) = DX(I+1)
       DY(I+2) = DX(I+2)
       DY(I+3) = DX(I+3)
       DY(I+4) = DX(I+4)
       DY(I+5) = DX(I+5)
       DY(I+6) = DX(I+6)
    end do
    RETURN
  end SUBROUTINE DCOPYN

  
#ifdef CUDA
  attributes(device) &
#endif
  SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
    INTEGER INCX,INCY,N
    DOUBLE PRECISION DX(:),DY(:)
! *  Purpose
! *  =======
! *
! *     copies a vector, x, to a vector, y.
! *     uses unrolled loops for increments equal to one.
! *     jack dongarra, linpack, 3/11/78.
! *     modified 12/3/93, array(1) declarations changed to array(*)
    INTEGER I,IX,IY,M,MP1

    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20

    ! code for unequal increments or equal increments
    ! not equal to 1

    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
       DY(IY) = DX(IX)
       IX = IX + INCX
       IY = IY + INCY
    end do
    RETURN

    ! code for both increments equal to 1

    ! clean-up loop

20  M = MOD(N,7)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
       DY(I) = DX(I)
    end do
    IF (N.LT.7) RETURN
40  MP1 = M + 1
    DO I = MP1,N,7
       DY(I) = DX(I)
       DY(I+1) = DX(I+1)
       DY(I+2) = DX(I+2)
       DY(I+3) = DX(I+3)
       DY(I+4) = DX(I+4)
       DY(I+5) = DX(I+5)
       DY(I+6) = DX(I+6)
    end do
    RETURN
  end SUBROUTINE DCOPY

#ifdef CUDA
  attributes(device) &
#endif  
  SUBROUTINE DAXPYN(N,DA,DX,INCX,DY,INCY)

    ! DAXPYN takes two vectors each of length N.
    ! This is not the case for the ordinary DAXPY.
  
    !$acc routine seq
    !     .. Scalar Arguments ..
    DOUBLE PRECISION DA
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(N),DY(N)
    !     ..
    !
    !  Purpose
    !   =======
    ! 
    !      constant times a vector plus a vector.
    !      uses unrolled loops for increments equal to one.
    !      jack dongarra, linpack, 3/11/78.
    !      modified 12/3/93, array(1) declarations changed to array(*)
    ! 
    ! 
    !      .. Local Scalars ..
    INTEGER I,IX,IY,M,MP1

    IF (N.LE.0) RETURN
    IF (DA.EQ.0.0d0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    ! 
    !         code for unequal increments or equal increments
    !           not equal to 1
    ! 
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
       DY(IY) = DY(IY) + DA*DX(IX)
       IX = IX + INCX
       IY = IY + INCY
    end do
    RETURN
    ! 
    !         code for both increments equal to 1
    ! 
    ! 
    !         clean-up loop
    ! 
20  M = MOD(N,4)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
       DY(I) = DY(I) + DA*DX(I)
    end do
    IF (N.LT.4) RETURN
40  MP1 = M + 1
    DO I = MP1,N,4
       DY(I) = DY(I) + DA*DX(I)
       DY(I+1) = DY(I+1) + DA*DX(I+1)
       DY(I+2) = DY(I+2) + DA*DX(I+2)
       DY(I+3) = DY(I+3) + DA*DX(I+3)
    end do
    RETURN
  END SUBROUTINE DAXPYN


#ifdef CUDA
  attributes(device) &
#endif  
  SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
    !$acc routine seq
    !     .. Scalar Arguments ..
    DOUBLE PRECISION DA
    INTEGER INCX,INCY,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DX(:),DY(:)
    !     ..
    !
    !  Purpose
    !   =======
    ! 
    !      constant times a vector plus a vector.
    !      uses unrolled loops for increments equal to one.
    !      jack dongarra, linpack, 3/11/78.
    !      modified 12/3/93, array(1) declarations changed to array(*)
    ! 
    ! 
    !      .. Local Scalars ..
    INTEGER I,IX,IY,M,MP1

    IF (N.LE.0) RETURN
    IF (DA.EQ.0.0d0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    ! 
    !         code for unequal increments or equal increments
    !           not equal to 1
    ! 
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
       DY(IY) = DY(IY) + DA*DX(IX)
       IX = IX + INCX
       IY = IY + INCY
    end do
    RETURN
    ! 
    !         code for both increments equal to 1
    ! 
    ! 
    !         clean-up loop
    ! 
20  M = MOD(N,4)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
       DY(I) = DY(I) + DA*DX(I)
    end do
    IF (N.LT.4) RETURN
40  MP1 = M + 1
    DO I = MP1,N,4
       DY(I) = DY(I) + DA*DX(I)
       DY(I+1) = DY(I+1) + DA*DX(I+1)
       DY(I+2) = DY(I+2) + DA*DX(I+2)
       DY(I+3) = DY(I+3) + DA*DX(I+3)
    end do
    RETURN
  END SUBROUTINE DAXPY
  
#ifdef CUDA
  attributes(device) &
#endif
  FUNCTION DDOT(N,DX,INCX,DY,INCY) result(dotval)
    DOUBLE PRECISION dotval
    !      .. Scalar Arguments ..
    INTEGER INCX,INCY,N
    !      ..
    !      .. Array Arguments ..
    DOUBLE PRECISION DX(:),DY(:)
    !      ..
    ! 
    !   Purpose
    !   =======
    ! 
    !      forms the dot product of two vectors.
    !      uses unrolled loops for increments equal to one.
    !      jack dongarra, linpack, 3/11/78.
    !      modified 12/3/93, array(1) declarations changed to array(*)
    ! 
    ! 
    !      .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY,M,MP1

    dotval = 0.0d0
    DTEMP = 0.0d0
    IF (N.LE.0) RETURN
    IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
    ! 
    !         code for unequal increments or equal increments
    !           not equal to 1
    ! 
    IX = 1
    IY = 1
    IF (INCX.LT.0) IX = (-N+1)*INCX + 1
    IF (INCY.LT.0) IY = (-N+1)*INCY + 1
    DO I = 1,N
       DTEMP = DTEMP + DX(IX)*DY(IY)
       IX = IX + INCX
       IY = IY + INCY
    end do
    dotval = DTEMP
    RETURN
    ! 
    !         code for both increments equal to 1
    ! 
    ! 
    !         clean-up loop
    ! 
20  M = MOD(N,5)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
       DTEMP = DTEMP + DX(I)*DY(I)
    end do
    IF (N.LT.5) GO TO 60
40  MP1 = M + 1
    DO I = MP1,N,5
       DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
    end do
60  dotval = DTEMP
    RETURN
  END FUNCTION DDOT

#ifdef CUDA
  attributes(device) &
#endif  
  SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    !$acc routine seq
    !      .. Scalar Arguments ..
    DOUBLE PRECISION ALPHA,BETA
    INTEGER K,LDA,LDB,LDC,M,N
    ! !      CHARACTER TRANSA,TRANSB
    INTEGER TRANSA,TRANSB
    !      ..
    !      .. Array Arguments ..
    DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
    !      ..
    ! 
    !   Purpose
    !   =======
    ! 
    !   DGEMM  performs one of the matrix-matrix operations
    ! 
    !      C := alpha*op( A )*op( B ) + beta*C,
    ! 
    !   where  op( X ) is one of
    ! 
    !      op( X ) = X   or   op( X ) = X**T,
    ! 
    !   alpha and beta are scalars, and A, B and C are matrices, with op( A )
    !   an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
    ! 
    !   Arguments
    !   ==========
    ! 
    !   TRANSA - CHARACTER*1.
    !            On entry, TRANSA specifies the form of op( A ) to be used in
    !            the matrix multiplication as follows:
    ! 
    !               TRANSA = 1 --> 'N' or 'n',  op( A ) = A.
    ! 
    !               TRANSA = 2 --> 'T' or 't',  op( A ) = A**T.
    ! 
    !               TRANSA = 3 --> 'C' or 'c',  op( A ) = A**T.
    ! 
    !            Unchanged on exit.
    ! 
    !   TRANSB - CHARACTER*1.
    !            On entry, TRANSB specifies the form of op( B ) to be used in
    !            the matrix multiplication as follows:
    ! 
    !               TRANSB = 1 --> 'N' or 'n',  op( B ) = B.
    !                              
    !               TRANSB = 2 --> 'T' or 't',  op( B ) = B**T.
    !                              
    !               TRANSB = 3 --> 'C' or 'c',  op( B ) = B**T.
    ! 
    !            Unchanged on exit.
    ! 
    !   M      - INTEGER.
    !            On entry,  M  specifies  the number  of rows  of the  matrix
    !            op( A )  and of the  matrix  C.  M  must  be at least  zero.
    !            Unchanged on exit.
    ! 
    !   N      - INTEGER.
    !            On entry,  N  specifies the number  of columns of the matrix
    !            op( B ) and the number of columns of the matrix C. N must be
    !            at least zero.
    !            Unchanged on exit.
    ! 
    !   K      - INTEGER.
    !            On entry,  K  specifies  the number of columns of the matrix
    !            op( A ) and the number of rows of the matrix op( B ). K must
    !            be at least  zero.
    !            Unchanged on exit.
    ! 
    !   ALPHA  - DOUBLE PRECISION.
    !            On entry, ALPHA specifies the scalar alpha.
    !            Unchanged on exit.
    ! 
    !   A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
    !            k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
    !            Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
    !            part of the array  A  must contain the matrix  A,  otherwise
    !            the leading  k by m  part of the array  A  must contain  the
    !            matrix A.
    !            Unchanged on exit.
    ! 
    !   LDA    - INTEGER.
    !            On entry, LDA specifies the first dimension of A as declared
    !            in the calling (sub) program. When  TRANSA = 'N' or 'n' then
    !            LDA must be at least  max( 1, m ), otherwise  LDA must be at
    !            least  max( 1, k ).
    !            Unchanged on exit.
    ! 
    !   B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
    !            n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
    !            Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
    !            part of the array  B  must contain the matrix  B,  otherwise
    !            the leading  n by k  part of the array  B  must contain  the
    !            matrix B.
    !            Unchanged on exit.
    ! 
    !   LDB    - INTEGER.
    !            On entry, LDB specifies the first dimension of B as declared
    !            in the calling (sub) program. When  TRANSB = 'N' or 'n' then
    !            LDB must be at least  max( 1, k ), otherwise  LDB must be at
    !            least  max( 1, n ).
    !            Unchanged on exit.
    ! 
    !   BETA   - DOUBLE PRECISION.
    !            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
    !            supplied as zero then C need not be set on input.
    !            Unchanged on exit.
    ! 
    !   C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
    !            Before entry, the leading  m by n  part of the array  C must
    !            contain the matrix  C,  except when  beta  is zero, in which
    !            case C need not be set on entry.
    !            On exit, the array  C  is overwritten by the  m by n  matrix
    !            ( alpha*op( A )*op( B ) + beta*C ).
    ! 
    !   LDC    - INTEGER.
    !            On entry, LDC specifies the first dimension of C as declared
    !            in  the  calling  (sub)  program.   LDC  must  be  at  least
    !            max( 1, m ).
    !            Unchanged on exit.
    ! 
    !   Further Details
    !   ===============
    ! 
    !   Level 3 Blas routine.
    ! 
    !   -- Written on 8-February-1989.
    !      Jack Dongarra, Argonne National Laboratory.
    !      Iain Duff, AERE Harwell.
    !      Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !      Sven Hammarling, Numerical Algorithms Group Ltd.
    ! 
    !   =====================================================================
    ! 
    !      .. External Functions ..
    ! LOGICAL LSAME
    ! EXTERNAL LSAME
    !      ..
    !      .. External Subroutines ..
    ! EXTERNAL XERBLA
    !      ..
    !      .. Local Scalars ..
    DOUBLE PRECISION TEMP
    INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
    LOGICAL NOTA,NOTB
    !      ..
    !      .. Parameters ..
    DOUBLE PRECISION ONE,ZERO
    PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
    !      ..
    ! 
    !      Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    !      transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    !      and  columns of  A  and the  number of  rows  of  B  respectively.
    ! 
    ! !      NOTA = LSAME(TRANSA,'N')
    NOTA = TRANSA == 1
    ! !      NOTB = LSAME(TRANSB,'N')
    NOTB = TRANSB == 1
    IF (NOTA) THEN
       NROWA = M
       NCOLA = K
    ELSE
       NROWA = K
       NCOLA = M
    END IF
    IF (NOTB) THEN
       NROWB = K
    ELSE
       NROWB = N
    END IF
    ! 
    !      Test the input parameters.
    ! 
    INFO = 0
    ! !      IF ((.NOT.NOTA) .AND. (.NOT.LSAME(TRANSA,'C')) .AND. &
    ! !         (.NOT.LSAME(TRANSA,'T'))) THEN
    IF ((.NOT.NOTA) .AND. (.NOT.TRANSA==3) .AND. &
         (.NOT.TRANSA==2)) THEN
       INFO = 1
       ! !      ELSE IF ((.NOT.NOTB) .AND. (.NOT.LSAME(TRANSB,'C')) .AND. &
       ! !              (.NOT.LSAME(TRANSB,'T'))) THEN
    ELSE IF ((.NOT.NOTB) .AND. (.NOT.TRANSB==3) .AND. &
         (.NOT.TRANSB==2)) THEN
       INFO = 2
    ELSE IF (M.LT.0) THEN
       INFO = 3
    ELSE IF (N.LT.0) THEN
       INFO = 4
    ELSE IF (K.LT.0) THEN
       INFO = 5
    ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
       INFO = 8
    ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
       INFO = 10
    ELSE IF (LDC.LT.MAX(1,M)) THEN
       INFO = 13
    END IF
    !IF (INFO.NE.0) THEN
    !    CALL XERBLA('DGEMM ',INFO)
    !    RETURN
    !END IF
    ! 
    !      Quick return if possible.
    ! 
    IF ((M.EQ.0) .OR. (N.EQ.0) .OR. &
         (((ALPHA.EQ.ZERO).OR. (K.EQ.0)).AND. (BETA.EQ.ONE))) RETURN
    ! 
    !      And if  alpha.eq.zero.
    ! 
    IF (ALPHA.EQ.ZERO) THEN
       IF (BETA.EQ.ZERO) THEN
          DO J = 1,N
             DO I = 1,M
                C(I,J) = ZERO
             end do
          end do
       ELSE
          DO J = 1,N
             DO I = 1,M
                C(I,J) = BETA*C(I,J)
             end do
          end do
       END IF
       RETURN
    END IF
    ! 
    !      Start the operations.
    ! 
    IF (NOTB) THEN
       IF (NOTA) THEN
          ! 
          !            Form  C := alpha*A*B + beta*C.
          ! 
          DO J = 1,N
             IF (BETA.EQ.ZERO) THEN
                DO I = 1,M
                   C(I,J) = ZERO
                end do
             ELSE IF (BETA.NE.ONE) THEN
                DO I = 1,M
                   C(I,J) = BETA*C(I,J)
                end do
             END IF
             DO L = 1,K
                IF (B(L,J).NE.ZERO) THEN
                   TEMP = ALPHA*B(L,J)
                   DO I = 1,M
                      C(I,J) = C(I,J) + TEMP*A(I,L)
                   end do
                END IF
             end do
          end do
       ELSE
          ! 
          !            Form  C := alpha*A**T*B + beta*C
          ! 
          DO J = 1,N
             DO I = 1,M
                TEMP = ZERO
                DO L = 1,K
                   TEMP = TEMP + A(L,I)*B(L,J)
                end do
                IF (BETA.EQ.ZERO) THEN
                   C(I,J) = ALPHA*TEMP
                ELSE
                   C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                END IF
             end do
          end do
       END IF
    ELSE
       IF (NOTA) THEN
          ! 
          !            Form  C := alpha*A*B**T + beta*C
          ! 
          DO J = 1,N
             IF (BETA.EQ.ZERO) THEN
                DO I = 1,M
                   C(I,J) = ZERO
                end do
             ELSE IF (BETA.NE.ONE) THEN
                DO I = 1,M
                   C(I,J) = BETA*C(I,J)
                end do
             END IF
             DO L = 1,K
                IF (B(J,L).NE.ZERO) THEN
                   TEMP = ALPHA*B(J,L)
                   DO I = 1,M
                      C(I,J) = C(I,J) + TEMP*A(I,L)
                   end do
                END IF
             end do
          end do
       ELSE
          ! 
          !            Form  C := alpha*A**T*B**T + beta*C
          ! 
          DO J = 1,N
             DO I = 1,M
                TEMP = ZERO
                DO L = 1,K
                   TEMP = TEMP + A(L,I)*B(J,L)
                end do
                IF (BETA.EQ.ZERO) THEN
                   C(I,J) = ALPHA*TEMP
                ELSE
                   C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                END IF
             end do
          end do
       END IF
    END IF
    ! 
    RETURN
    ! 
    !      End of DGEMM .
    ! 
  end SUBROUTINE DGEMM

#ifdef CUDA
  attributes(device) &
#endif    
  SUBROUTINE DSCALN(N,DA,DX,INCX)

    ! DSCALN takes one vector of length N.
    ! This is not the case for the ordinary DSCALN.

    !$acc routine seq
    !      .. Scalar Arguments ..
    DOUBLE PRECISION DA
    INTEGER INCX,N
    !      ..
    !      .. Array Arguments ..
    DOUBLE PRECISION DX(N)
    !      ..
    ! 
    !   Purpose
    !   =======
    ! *
    !      scales a vector by a constant.
    !      uses unrolled loops for increment equal to one.
    !      jack dongarra, linpack, 3/11/78.
    !      modified 3/93 to return if incx .le. 0.
    !      modified 12/3/93, array(1) declarations changed to array(*)
    ! 
    ! 
    !      .. Local Scalars ..
    INTEGER I,M,MP1,NINCX
    !      ..
    !      .. Intrinsic Functions ..
    INTRINSIC MOD
    !      ..
    IF (N.LE.0 .OR. INCX.LE.0) RETURN
    IF (INCX.EQ.1) GO TO 20
    ! 
    !         code for increment not equal to 1
    ! 
    NINCX = N*INCX
    DO I = 1,NINCX,INCX
       DX(I) = DA*DX(I)
    end do
    RETURN
    ! 
    !         code for increment equal to 1
    ! 
    ! 
    !         clean-up loop
    ! 
20  M = MOD(N,5)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
       DX(I) = DA*DX(I)
    end do
    IF (N.LT.5) RETURN
40  MP1 = M + 1
    DO I = MP1,N,5
       DX(I) = DA*DX(I)
       DX(I+1) = DA*DX(I+1)
       DX(I+2) = DA*DX(I+2)
       DX(I+3) = DA*DX(I+3)
       DX(I+4) = DA*DX(I+4)
    end do
    RETURN
  END SUBROUTINE DSCALN


#ifdef CUDA
  attributes(device) &
#endif    
  SUBROUTINE DSCAL(N,DA,DX,INCX)
    !$acc routine seq
    !      .. Scalar Arguments ..
    DOUBLE PRECISION DA
    INTEGER INCX,N
    !      ..
    !      .. Array Arguments ..
    DOUBLE PRECISION DX(:)
    !      ..
    ! 
    !   Purpose
    !   =======
    ! *
    !      scales a vector by a constant.
    !      uses unrolled loops for increment equal to one.
    !      jack dongarra, linpack, 3/11/78.
    !      modified 3/93 to return if incx .le. 0.
    !      modified 12/3/93, array(1) declarations changed to array(*)
    ! 
    ! 
    !      .. Local Scalars ..
    INTEGER I,M,MP1,NINCX
    !      ..
    !      .. Intrinsic Functions ..
    INTRINSIC MOD
    !      ..
    IF (N.LE.0 .OR. INCX.LE.0) RETURN
    IF (INCX.EQ.1) GO TO 20
    ! 
    !         code for increment not equal to 1
    ! 
    NINCX = N*INCX
    DO I = 1,NINCX,INCX
       DX(I) = DA*DX(I)
    end do
    RETURN
    ! 
    !         code for increment equal to 1
    ! 
    ! 
    !         clean-up loop
    ! 
20  M = MOD(N,5)
    IF (M.EQ.0) GO TO 40
    DO I = 1,M
       DX(I) = DA*DX(I)
    end do
    IF (N.LT.5) RETURN
40  MP1 = M + 1
    DO I = MP1,N,5
       DX(I) = DA*DX(I)
       DX(I+1) = DA*DX(I+1)
       DX(I+2) = DA*DX(I+2)
       DX(I+3) = DA*DX(I+3)
       DX(I+4) = DA*DX(I+4)
    end do
    RETURN
  END SUBROUTINE DSCAL

#ifdef CUDA
  attributes(device) &
#endif    
  INTEGER FUNCTION IDAMAX(N,DX,INCX)
    !$acc routine seq
    !      .. Scalar Arguments ..
    INTEGER INCX,N
    !      ..
    !      .. Array Arguments ..
    DOUBLE PRECISION DX(:)
    !      ..
    ! 
    !   Purpose
    !   =======
    ! 
    !      finds the index of element having max. absolute value.
    !      jack dongarra, linpack, 3/11/78.
    !      modified 3/93 to return if incx .le. 0.
    !      modified 12/3/93, array(1) declarations changed to array(*)
    ! 
    ! 
    !      .. Local Scalars ..
    DOUBLE PRECISION DMAX
    INTEGER I,IX
    !      ..
    !      .. Intrinsic Functions ..
    INTRINSIC DABS
    !      ..
    IDAMAX = 0
    IF (N.LT.1 .OR. INCX.LE.0) RETURN
    IDAMAX = 1
    IF (N.EQ.1) RETURN
    IF (INCX.EQ.1) GO TO 20
    ! 
    !         code for increment not equal to 1
    ! 
    IX = 1
    DMAX = DABS(DX(1))
    IX = IX + INCX
    DO I = 2,N
       IF (DABS(DX(IX)).LE.DMAX) GO TO 5
       IDAMAX = I
       DMAX = DABS(DX(IX))
5      IX = IX + INCX
    end do
    RETURN
    ! 
    !         code for increment equal to 1
    ! 
20  DMAX = DABS(DX(1))
    DO I = 2,N
       IF (DABS(DX(I)).LE.DMAX) cycle
       IDAMAX = I
       DMAX = DABS(DX(I))
    end do
    RETURN
  END FUNCTION IDAMAX
  
end module blas_module
