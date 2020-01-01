module linpack_module

  use cuvode_parameters_module, only: VODE_NEQS

  use amrex_fort_module, only : rt => amrex_real
  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dgesl(a, ipvt, b)

    implicit none

    integer, parameter :: lda = VODE_NEQS
    integer, parameter :: n = VODE_NEQS
    integer, parameter :: job = 0

    integer, intent(in) :: ipvt(n)
    real(rt)        , intent(in) :: a(lda,n)
    real(rt)        , intent(inout) :: b(n)

    ! 
    !      dgesl solves the real(rt) system
    !      a * x = b  or  trans(a) * x = b
    !      using the factors computed by dgeco or dgefa.
    !
    !      on entry
    ! 
    !         a       real(rt)        (lda, n)
    !                 the output from dgeco or dgefa.
    !
    !         lda     integer
    !                 the leading dimension of the array  a .
    !
    !         n       integer
    !                 the order of the matrix  a .
    !
    !         ipvt    integer(n)
    !                 the pivot array from dgeco or dgefa.
    ! 
    !         b       real(rt)        (n)
    !                 the right hand side array.
    !
    !         job     integer
    !                 = 0         to solve  a*x = b ,
    !                 = nonzero   to solve  trans(a)*x = b  where
    !                             trans(a)  is the transpose.
    !
    !      on return
    !
    !         b       the solution array  x .
    !
    !      error condition
    !
    !         a division by zero will occur if the input factor contains a
    !         zero on the diagonal.  technically this indicates singularity
    !         but it is often caused by improper arguments or improper
    !         setting of lda .  it will not occur if the subroutines are
    !         called correctly and if dgeco has set rcond .gt. 0.0
    !         or dgefa has set info .eq. 0 .
    !
    !      to compute  inverse(a) * c  where  c  is a matrix
    !      with  p  columns
    !            call dgeco(a,lda,n,ipvt,rcond,z)
    !            if (rcond is too small) go to ...
    !            do 10 j = 1, p
    !               call dgesl(a,lda,n,ipvt,c(1,j),0)
    !         10 continue
    !
    !      linpack. this version dated 08/14/78 .
    !      cleve moler, university of new mexico, argonne national lab.
    !
    !      internal variables
    ! 
    real(rt)         t

    integer k,kb,l,nm1
    !$gpu
    !
    nm1 = n - 1
    if (job .ne. 0) goto 50
    !
    !         job = 0 , solve  a * x = b
    !         first solve  l*y = b
    !
    if (nm1 .lt. 1) goto 30
    do k = 1, nm1
       l = ipvt(k)
       t = b(l)
       if (l .eq. k) goto 10
       b(l) = b(k)
       b(k) = t
10     continue
       b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)
    enddo
30  continue
    !
    !         now solve  u*x = y
    !
    do kb = 1, n
       k = n + 1 - kb
       b(k) = b(k)/a(k,k)
       t = -b(k)
       b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
    enddo
    go to 100
50  continue
    !
    !         job = nonzero, solve  trans(a) * x = b
    !         first solve  trans(u)*y = b
    !
    do k = 1, n
       t = sum(a(1:k-1,k) * b(1:k-1))
       b(k) = (b(k) - t)/a(k,k)
    enddo
    !
    !         now solve trans(l)*x = y
    !
    if (nm1 .lt. 1) goto 90
    do kb = 1, nm1
       k = n - kb
       b(k) = b(k) + sum(a(k+1:n,k) * b(k+1:n))
       l = ipvt(k)
       if (l .eq. k) go to 70
       t = b(l)
       b(l) = b(k)
       b(k) = t
70     continue
    enddo
90  continue
100 continue
  end subroutine dgesl

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine DGBFA (ABD, LDA, ML, MU, IPVT, INFO)

    ! ***BEGIN PROLOGUE  DGBFA
    ! ***PURPOSE  Factor a band matrix using Gaussian elimination.
    ! ***CATEGORY  D2A2
    ! ***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
    ! ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
    ! ***AUTHOR  Moler, C. B., (U. of New Mexico)
    ! ***DESCRIPTION
    !
    !      DGBFA factors a double precision band matrix by elimination.
    !
    !      DGBFA is usually called by DGBCO, but it can be called
    !      directly with a saving in time if  RCOND  is not needed.
    !
    !      On Entry
    !
    !         ABD     DOUBLE PRECISION(LDA, N)
    !                 contains the matrix in band storage.  The columns
    !                 of the matrix are stored in the columns of  ABD  and
    !                 the diagonals of the matrix are stored in rows
    !                 ML+1 through 2*ML+MU+1 of  ABD .
    !                 See the comments below for details.
    !
    !         LDA     INTEGER
    !                 the leading dimension of the array  ABD .
    !                 LDA must be .GE. 2*ML + MU + 1 .
    !
    !         N       INTEGER
    !                 the order of the original matrix.
    !
    !         ML      INTEGER
    !                 number of diagonals below the main diagonal.
    !                 0 .LE. ML .LT.  N .
    !
    !         MU      INTEGER
    !                 number of diagonals above the main diagonal.
    !                 0 .LE. MU .LT.  N .
    !                 More efficient if  ML .LE. MU .
    !      On Return
    !
    !         ABD     an upper triangular matrix in band storage and
    !                 the multipliers which were used to obtain it.
    !                 The factorization can be written  A = L*U  where
    !                 L  is a product of permutation and unit lower
    !                 triangular matrices and  U  is upper triangular.
    !
    !         IPVT    INTEGER(N)
    !                 an integer array of pivot indices.
    !
    !         INFO    INTEGER
    !                 = 0  normal value.
    !                 = K  if  U(K,K) .EQ. 0.0 .  This is not an error
    !                      condition for this subroutine, but it does
    !                      indicate that DGBSL will divide by zero if
    !                      called.  Use  RCOND  in DGBCO for a reliable
    !                      indication of singularity.
    !
    !      Band Storage
    !
    !            If  A  is a band matrix, the following program segment
    !            will set up the input.
    !
    !                    ML = (band width below the diagonal)
    !                    MU = (band width above the diagonal)
    !                    M = ML + MU + 1
    !                    DO 20 J = 1, N
    !                       I1 = MAX(1, J-MU)
    !                       I2 = MIN(N, J+ML)
    !                       DO 10 I = I1, I2
    !                          K = I - J + M
    !                          ABD(K,J) = A(I,J)
    !                 10    CONTINUE
    !                 20 CONTINUE
    !
    !            This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
    !            In addition, the first  ML  rows in  ABD  are used for
    !            elements generated during the triangularization.
    !            The total number of rows needed in  ABD  is  2*ML+MU+1 .
    !            The  ML+MU by ML+MU  upper left triangle and the
    !            ML by ML  lower right triangle are not referenced.
    !
    ! ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
    !                  Stewart, LINPACK Users' Guide, SIAM, 1979.
    ! ***REVISION HISTORY  (YYMMDD)
    !    780814  DATE WRITTEN
    !    890531  Changed all specific intrinsics to generic.  (WRB)
    !    890831  Modified array declarations.  (WRB)
    !    890831  REVISION DATE from Version 3.2
    !    891214  Prologue converted to Version 4.0 format.  (BAB)
    !    900326  Removed duplicate information from DESCRIPTION section.
    !            (WRB)
    !    920501  Reformatted the REFERENCES section.  (WRB)
  ! ***END PROLOGUE  DGBFA
    integer, parameter :: N = VODE_NEQS
    INTEGER LDA,ML,MU,IPVT(:),INFO
    real(rt)         ABD(LDA, N), dABD(LDA)
    ! 
    real(rt)         T

    INTEGER I,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
    !$gpu
    !
    ! ***FIRST EXECUTABLE STATEMENT  DGBFA
    M = ML + MU + 1
    INFO = 0
    !
    !      ZERO INITIAL FILL-IN COLUMNS
    !
    J0 = MU + 2
    J1 = MIN(N,M) - 1
    IF (J1 .LT. J0) GO TO 30
    DO JZ = J0, J1
       I0 = M + 1 - JZ
       DO I = I0, ML
          ABD(I,JZ) = 0.0e0_rt
       end do
    end do
30  CONTINUE
    JZ = J1
    JU = 0
    !
    !      GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
    !
    NM1 = N - 1
    IF (NM1 .LT. 1) GO TO 130
    DO K = 1, NM1
       KP1 = K + 1
       !
       !         ZERO NEXT FILL-IN COLUMN
       !
       JZ = JZ + 1
       IF (JZ .GT. N) GO TO 50
       IF (ML .LT. 1) GO TO 50
       DO I = 1, ML
          ABD(I,JZ) = 0.0e0_rt
       end do
50     CONTINUE
       !
       !         FIND L = PIVOT INDEX
       !
       LM = MIN(ML,N-K)
       L = idamax(LM+1,ABD(M:M+LM,K)) + M - 1
       IPVT(K) = L + K - M
       !
       !         ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
       ! 
       IF (ABD(L,K) .EQ. 0.0e0_rt) GO TO 100
       ! 
       !            INTERCHANGE IF NECESSARY
       !
       IF (L .EQ. M) GO TO 60
       T = ABD(L,K)
       ABD(L,K) = ABD(M,K)
       ABD(M,K) = T
60     CONTINUE
       !
       !            COMPUTE MULTIPLIERS
       ! 
       T = -1.0e0_rt/ABD(M,K)
       ABD(M+1:M+LM,K) = ABD(M+1:M+LM,K) * T
       !
       !            ROW ELIMINATION WITH COLUMN INDEXING
       !
       JU = MIN(MAX(JU,MU+IPVT(K)),N)
       MM = M
       IF (JU .LT. KP1) GO TO 90
       DO J = KP1, JU
          L = L - 1
          MM = MM - 1
          T = ABD(L,J)
          IF (L .EQ. MM) GO TO 70
          ABD(L,J) = ABD(MM,J)
          ABD(MM,J) = T
70        CONTINUE
          dABD(M+1:M+LM) = T * ABD(M+1:M+LM,K)
          ABD(MM+1:MM+LM,J) = ABD(MM+1:MM+LM,J) + dABD(M+1:M+LM)
       end do
90     CONTINUE
       GO TO 110
100    CONTINUE
       INFO = K
110    CONTINUE
    end do
130 CONTINUE
    IPVT(N) = N
    IF (ABD(M,N) .EQ. 0.0e0_rt) INFO = N
    RETURN
  END SUBROUTINE DGBFA

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  SUBROUTINE DGBSL (ABD, LDA, ML, MU, IPVT, B, JOB)

    ! ***BEGIN PROLOGUE  DGBSL
    ! ***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
    !             the factors computed by DGBCO or DGBFA.
    ! ***CATEGORY  D2A2
    ! ***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
    ! ***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
    ! ***AUTHOR  Moler, C. B., (U. of New Mexico)
    ! ***DESCRIPTION
    !
    !      DGBSL solves the double precision band system
    !      A * X = B  or  TRANS(A) * X = B
    !      using the factors computed by DGBCO or DGBFA.
    !
    !      On Entry
    !
    !         ABD     DOUBLE PRECISION(LDA, N)
    !                 the output from DGBCO or DGBFA.
    !
    !         LDA     INTEGER
    !                 the leading dimension of the array  ABD .
    !
    !         N       INTEGER
    !                 the order of the original matrix.
    !
    !         ML      INTEGER
    !                 number of diagonals below the main diagonal.
    !
    !         MU      INTEGER
    !                 number of diagonals above the main diagonal.
    !
    !         IPVT    INTEGER(N)
    !                 the pivot array from DGBCO or DGBFA.
    !
    !         B       DOUBLE PRECISION(N)
    !                 the right hand side array.
    !
    !         JOB     INTEGER
    !                 = 0         to solve  A*X = B ,
    !                 = nonzero   to solve  TRANS(A)*X = B , where
    !                             TRANS(A)  is the transpose.
    !
    !      On Return
    !
    !         B       the solution array  X .
    !
    !      Error Condition
    !
    !         A division by zero will occur if the input factor contains a
    !         zero on the diagonal.  Technically this indicates singularity
    !         but it is often caused by improper arguments or improper
    !         setting of LDA .  It will not occur if the subroutines are
    !         called correctly and if DGBCO has set RCOND .GT. 0.0
    !         or DGBFA has set INFO .EQ. 0 .
    !
    !      To compute  INVERSE(A) * C  where  C  is a matrix
    !      with  P  columns
    !            CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
    !            IF (RCOND is too small) GO TO ...
    !            DO 10 J = 1, P
    !               CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
    !         10 CONTINUE
    !
    ! ***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
    !                  Stewart, LINPACK Users' Guide, SIAM, 1979.
    ! ***REVISION HISTORY  (YYMMDD)
    !    780814  DATE WRITTEN
    !    890531  Changed all specific intrinsics to generic.  (WRB)
    !    890831  Modified array declarations.  (WRB)
    !    890831  REVISION DATE from Version 3.2
    !    891214  Prologue converted to Version 4.0 format.  (BAB)
    !    900326  Removed duplicate information from DESCRIPTION section.
    !            (WRB)
    !    920501  Reformatted the REFERENCES section.  (WRB)
    ! ***END PROLOGUE  DGBSL
    integer, parameter :: N = VODE_NEQS
    INTEGER LDA,ML,MU,IPVT(:),JOB
    real(rt)         ABD(LDA,N),B(:)
    ! 
    real(rt)         T
    INTEGER K,KB,L,LA,LB,LM,M,NM1
    !$gpu
    ! ***FIRST EXECUTABLE STATEMENT  DGBSL
    M = MU + ML + 1
    NM1 = N - 1
    IF (JOB .NE. 0) GO TO 50
    !
    !         JOB = 0 , SOLVE  A * X = B
    !         FIRST SOLVE L*Y = B
    !
    IF (ML .EQ. 0) GO TO 30
    IF (NM1 .LT. 1) GO TO 30
    DO K = 1, NM1
       LM = MIN(ML,N-K)
       L = IPVT(K)
       T = B(L)
       IF (L .EQ. K) GO TO 10
       B(L) = B(K)
       B(K) = T
10     CONTINUE
       B(K+1:K+LM) = B(K+1:K+LM) + T * ABD(M+1:M+LM,K)
    end do
30  CONTINUE
    !
    !         NOW SOLVE  U*X = Y
    !
    DO KB = 1, N
       K = N + 1 - KB
       B(K) = B(K)/ABD(M,K)
       LM = MIN(K,M) - 1
       LA = M - LM
       LB = K - LM
       T = -B(K)
       B(LB:LB+LM-1) = B(LB:LB+LM-1) + T * ABD(LA:LA+LM-1,K)
    end do
    GO TO 100
50  CONTINUE
    !
    !         JOB = NONZERO, SOLVE  TRANS(A) * X = B
    !         FIRST SOLVE  TRANS(U)*Y = B
    !
    DO K = 1, N
       LM = MIN(K,M) - 1
       LA = M - LM
       LB = K - LM
       T = sum(ABD(LA:LA + LM - 1,K) * B(LB:LB + LM - 1))
       B(K) = (B(K) - T)/ABD(M,K)
    end do
    !
    !         NOW SOLVE TRANS(L)*X = Y
    !
    IF (ML .EQ. 0) GO TO 90
    IF (NM1 .LT. 1) GO TO 90
    DO KB = 1, NM1
       K = N - KB
       LM = MIN(ML,N-K)
       B(K) = B(K) + sum(ABD(M+1:M+LM,K) * B(K+1:K+LM))
       L = IPVT(K)
       IF (L .EQ. K) GO TO 70
       T = B(L)
       B(L) = B(K)
       B(K) = T
70     CONTINUE
    end do
90  CONTINUE
100 CONTINUE
    RETURN
  END SUBROUTINE DGBSL

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dgefa (a, ipvt, info)

    integer, parameter :: lda = VODE_NEQS
    integer, parameter :: n = VODE_NEQS
    integer ipvt(n), info
    real(rt)         a(lda, n)
    ! 
    !      dgefa factors a real(rt) matrix by gaussian elimination.
    ! 
    !      dgefa is usually called by dgeco, but it can be called
    !      directly with a saving in time if  rcond  is not needed.
    !      (time for dgeco) = (1 + 9/n)*(time for dgefa) .
    !
    !      on entry
    ! 
    !         a       real(rt)        (lda, n)
    !                 the matrix to be factored.
    !
    !         lda     integer
    !                 the leading dimension of the array  a .
    !
    !         n       integer
    !                 the order of the matrix  a .
    !
    !      on return
    !
    !         a       an upper triangular matrix and the multipliers
    !                 which were used to obtain it.
    !                 the factorization can be written  a = l*u  where
    !                 l  is a product of permutation and unit lower
    !                 triangular matrices and  u  is upper triangular.
    !
    !         ipvt    integer(n)
    !                 an integer array of pivot indices.
    !
    !         info    integer
    !                 = 0  normal value.
    !                 = k  if  u(k,k) .eq. 0.0 .  this is not an error
    !                      condition for this subroutine, but it does
    !                      indicate that dgesl or dgedi will divide by zero
    !                      if called.  use  rcond  in dgeco for a reliable
    !                      indication of singularity.
    !
    !      linpack. this version dated 08/14/78 .
    !      cleve moler, university of new mexico, argonne national lab.
    !
    !      subroutines and functions
    !
    !      blas vdaxpy
    !
    !      internal variables
    ! 
    real(rt)         t
    integer j,k,kp1,l,nm1
    !$gpu
    !
    !
    !      gaussian elimination with partial pivoting
    !

    info = 0
    nm1 = n - 1

    if (nm1 .lt. 1) goto 70
    do k = 1, nm1
       kp1 = k + 1
       !
       !         find l = pivot index
       !
       l = idamax(n-k+1,a(k:n,k)) + k - 1
       ipvt(k) = l
       !
       !         zero pivot implies this column already triangularized
       ! 
       if (a(l,k) .eq. 0.0e0_rt) goto 40
       ! 
       !            interchange if necessary
       !
       if (l .eq. k) goto 10
       t = a(l,k)
       a(l,k) = a(k,k)
       a(k,k) = t
10     continue
       !
       !            compute multipliers
       ! 
       t = -1.0e0_rt/a(k,k)
       a(k+1:n,k) = a(k+1:n,k) * t
       !
       !            row elimination with column indexing
       !
       do j = kp1, n
          t = a(l,j)
          if (l .eq. k) goto 20
          a(l,j) = a(k,j)
          a(k,j) = t
20        continue
          a(k+1:n,j) = a(k+1:n,j) + t * a(k+1:n,k)
       enddo
       goto 50
40     continue
       info = k
50     continue
    enddo
70  continue
    ipvt(n) = n
    if (a(n,n) .eq. 0.0e0_rt) info = n
  end subroutine dgefa

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  function idamax(N, x) result(index)

    implicit none

    integer, intent(in) :: N
    real(rt)        , intent(in) :: x(N)

    real(rt)         :: dmax
    integer :: i, index

    !$gpu

    index = 1

    dmax = abs(x(1))
    DO i = 2, N
       IF (abs(X(i)) .le. dmax) cycle
       index = i
       dmax = abs(x(i))
    end do
  end function idamax

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

    use amrex_constants_module, only : ZERO, ONE

    implicit none

    ! scalar arguments
    double precision, intent(in) :: alpha, beta
    integer, intent(in) :: k, lda, ldb, ldc, m, n

    ! character transa,transb
    integer, intent(in) :: transa, transb

    !  array arguments
    double precision :: a(lda,*), b(ldb,*), c(ldc,*)

    !  Purpose
    !  =======
    !
    !  DGEMM  performs one of the matrix-matrix operations
    !
    !     C := alpha*op( A )*op( B ) + beta*C,
    !
    !  where  op( X ) is one of
    !
    !     op( X ) = X   or   op( X ) = X**T,
    !
    !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
    !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
    !
    !  Arguments
    !  ==========
    !
    !  TRANSA - CHARACTER*1.
    !           On entry, TRANSA specifies the form of op( A ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSA = 1 --> 'N' or 'n',  op( A ) = A.
    !
    !              TRANSA = 2 --> 'T' or 't',  op( A ) = A**T.
    !
    !              TRANSA = 3 --> 'C' or 'c',  op( A ) = A**T.
    !
    !           Unchanged on exit.
    !
    !  TRANSB - CHARACTER*1.
    !           On entry, TRANSB specifies the form of op( B ) to be used in
    !           the matrix multiplication as follows:
    !
    !              TRANSB = 1 --> 'N' or 'n',  op( B ) = B.
    !
    !              TRANSB = 2 --> 'T' or 't',  op( B ) = B**T.
    !
    !              TRANSB = 3 --> 'C' or 'c',  op( B ) = B**T.
    !
    !           Unchanged on exit.
    !
    !  M      - INTEGER.
    !           On entry,  M  specifies  the number  of rows  of the  matrix
    !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry,  N  specifies the number  of columns of the matrix
    !           op( B ) and the number of columns of the matrix C. N must be
    !           at least zero.
    !           Unchanged on exit.
    !
    !  K      - INTEGER.
    !           On entry,  K  specifies  the number of columns of the matrix
    !           op( A ) and the number of rows of the matrix op( B ). K must
    !           be at least  zero.
    !           Unchanged on exit.
    !
    !  ALPHA  - DOUBLE PRECISION.
    !           On entry, ALPHA specifies the scalar alpha.
    !           Unchanged on exit.
    !
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
    !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
    !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
    !           part of the array  A  must contain the matrix  A,  otherwise
    !           the leading  k by m  part of the array  A  must contain  the
    !           matrix A.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
    !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
    !           least  max( 1, k ).
    !           Unchanged on exit.
    !
    !  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
    !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
    !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
    !           part of the array  B  must contain the matrix  B,  otherwise
    !           the leading  n by k  part of the array  B  must contain  the
    !           matrix B.
    !           Unchanged on exit.
    !
    !  LDB    - INTEGER.
    !           On entry, LDB specifies the first dimension of B as declared
    !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
    !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
    !           least  max( 1, n ).
    !           Unchanged on exit.
    !
    !  BETA   - DOUBLE PRECISION.
    !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
    !           supplied as zero then C need not be set on input.
    !           Unchanged on exit.
    !
    !  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
    !           Before entry, the leading  m by n  part of the array  C must
    !           contain the matrix  C,  except when  beta  is zero, in which
    !           case C need not be set on entry.
    !           On exit, the array  C  is overwritten by the  m by n  matrix
    !           ( alpha*op( A )*op( B ) + beta*C ).
    !
    !  LDC    - INTEGER.
    !           On entry, LDC specifies the first dimension of C as declared
    !           in  the  calling  (sub)  program.   LDC  must  be  at  least
    !           max( 1, m ).
    !           Unchanged on exit.
    !
    !  Further Details
    !  ===============
    !
    !  Level 3 Blas routine.
    !
    !  -- Written on 8-February-1989.
    !     Jack Dongarra, Argonne National Laboratory.
    !     Iain Duff, AERE Harwell.
    !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !     Sven Hammarling, Numerical Algorithms Group Ltd.
    !
    !  =====================================================================
    !


    ! Local Scalars ..
    double precision temp
    integer i, info, j, l, ncola, nrowa, nrowb
    logical nota, notb

    !$gpu

    !
    !
    ! Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
    ! transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
    ! and  columns of  A  and the  number of  rows  of  B  respectively.

    nota = transa == 1
    notb = transb == 1
    if (nota) then
       nrowa = m
       ncola = k
    else
       nrowa = k
       ncola = m
    end if

    if (notb) then
       nrowb = k
    else
       nrowb = n
    end if

    ! Test the input parameters.
    info = 0
    if ((.not.nota) .and. (.not.transa==3) .and. (.not.transa==2)) then
       info = 1
    else if ((.not.notb) .and. (.not.transb==3) .and. (.not.transb==2)) then
       info = 2
    else if (m.lt.0) then
       info = 3
    else if (n.lt.0) then
       info = 4
    else if (k.lt.0) then
       info = 5
    else if (lda.lt.max(1,nrowa)) then
       info = 8
    else if (ldb.lt.max(1,nrowb)) then
       info = 10
    else if (ldc.lt.max(1,m)) then
       info = 13
    end if

    !  Quick return if possible.

    if ((m .eq. 0) .or. (n .eq. 0) .or. &
         (((alpha .eq. ZERO) .or. (k .eq. 0)) .and. (beta .eq. ONE))) return

    ! And if  alpha.eq.zero.

    if (alpha .eq. zero) then
       if (beta .eq. zero) then
          do j = 1,n
             do i = 1,m
                c(i,j) = zero
             end do
          end do

       else
          do j = 1,n
             do i = 1,m
                c(i,j) = beta*c(i,j)
             end do
          end do
       end if
       return
    end if

    ! Start the operations.

    if (notb) then
       if (nota) then

          ! Form  C := alpha*A*B + beta*C.

          do j = 1,n
             if (beta.eq.zero) then
                do i = 1,m
                   c(i,j) = zero
                end do
             else if (beta.ne.one) then
                do i = 1,m
                   c(i,j) = beta*c(i,j)
                end do
             end if
             do l = 1,k
                if (b(l,j).ne.zero) then
                   temp = alpha*b(l,j)
                   do i = 1,m
                      c(i,j) = c(i,j) + temp*a(i,l)
                   end do
                end if
             end do
          end do

       else

          ! Form  C := alpha*A**T*B + beta*C

          do j = 1,n
             do i = 1,m
                temp = zero
                do l = 1,k
                   temp = temp + a(l,i)*b(l,j)
                end do
                if (beta.eq.zero) then
                   c(i,j) = alpha*temp
                else
                   c(i,j) = alpha*temp + beta*c(i,j)
                end if
             end do
          end do
       end if

    else

       if (nota) then

          ! Form  C := alpha*A*B**T + beta*C

              do j = 1,n
                 if (beta.eq.zero) then
                    do i = 1,m
                       c(i,j) = zero
                    end do
                 else if (beta.ne.one) then
                    do i = 1,m
                       c(i,j) = beta*c(i,j)
                    end do
                 end if
                 do l = 1,k
                    if (b(j,l).ne.zero) then
                       temp = alpha*b(j,l)
                       do i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
                       end do
                    end if
                 end do
              end do

          else

             ! form  c := alpha*a**t*b**t + beta*c

             do j = 1,n
                do i = 1,m
                   temp = zero
                   do l = 1,k
                      temp = temp + a(l,i)*b(j,l)
                   end do
                   if (beta.eq.zero) then
                      c(i,j) = alpha*temp
                   else
                      c(i,j) = alpha*temp + beta*c(i,j)
                   end if
                end do
             end do
          end if
       end if

     end subroutine dgemm

end module linpack_module
