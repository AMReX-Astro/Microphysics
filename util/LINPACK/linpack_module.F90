module linpack_module

  use blas_module
  
  implicit none
  
contains
  
#ifdef CUDA
  attributes(device) &
#endif
  subroutine dgesl (a,lda,n,ipvt,b,job)

    !$acc routine seq
    !$acc routine(daxpy) seq
    !$acc routine(vddot) seq

    integer lda,n,ipvt(:),job
    double precision a(lda,n),b(:)
    ! 
    !      dgesl solves the double precision system
    !      a * x = b  or  trans(a) * x = b
    !      using the factors computed by dgeco or dgefa.
    ! 
    !      on entry
    ! 
    !         a       double precision(lda, n)
    !                 the output from dgeco or dgefa.
    ! 
    !         lda     integer
    !                 the leading dimension of the array  a .
    ! 
    !         n       integer
    !                 the order of the matrix  a .
    ! 
    !         ipvt    integer(n)
    !                 the pivot vector from dgeco or dgefa.
    ! 
    !         b       double precision(n)
    !                 the right hand side vector.
    ! 
    !         job     integer
    !                 = 0         to solve  a*x = b ,
    !                 = nonzero   to solve  trans(a)*x = b  where
    !                             trans(a)  is the transpose.
    ! 
    !      on return
    ! 
    !         b       the solution vector  x .
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
    !      subroutines and functions
    ! 
    !      blas vdaxpy,vddot
    ! 
    !      internal variables
    ! 
    double precision t
    integer k,kb,l,nm1
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
       call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
    enddo
30  continue
    ! 
    !         now solve  u*x = y
    ! 
    do kb = 1, n
       k = n + 1 - kb
       b(k) = b(k)/a(k,k)
       t = -b(k)
       call daxpy(k-1,t,a(1,k),1,b(1),1)
    enddo
    go to 100
50  continue
    ! 
    !         job = nonzero, solve  trans(a) * x = b
    !         first solve  trans(u)*y = b
    ! 
    do k = 1, n
       t = vddot(k-1,a(1,k),1,b(1),1)
       b(k) = (b(k) - t)/a(k,k)
    enddo
    ! 
    !         now solve trans(l)*x = y
    ! 
    if (nm1 .lt. 1) goto 90
    do kb = 1, nm1
       k = n - kb
       b(k) = b(k) + vddot(n-k,a(k+1,k),1,b(k+1),1)
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

#ifdef CUDA
  attributes(device) &
#endif  
  SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
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
    !                 an integer vector of pivot indices.
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
    ! ***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
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
    INTEGER LDA,N,ML,MU,IPVT(:),INFO
    DOUBLE PRECISION ABD(:,:)
    ! 
    DOUBLE PRECISION T
    INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
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
          ABD(I,JZ) = 0.0D0
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
          ABD(I,JZ) = 0.0D0
       end do
50     CONTINUE
       ! 
       !         FIND L = PIVOT INDEX
       ! 
       LM = MIN(ML,N-K)
       L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
       IPVT(K) = L + K - M
       ! 
       !         ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
       ! 
       IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
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
       T = -1.0D0/ABD(M,K)
       CALL DSCAL(LM,T,ABD(M+1,K),1)
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
          CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
       end do
90     CONTINUE
       GO TO 110
100    CONTINUE
       INFO = K
110    CONTINUE
    end do
130 CONTINUE
    IPVT(N) = N
    IF (ABD(M,N) .EQ. 0.0D0) INFO = N
    RETURN
  END SUBROUTINE DGBFA

#ifdef CUDA
  attributes(device) &
#endif  
  SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
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
    !                 the pivot vector from DGBCO or DGBFA.
    ! 
    !         B       DOUBLE PRECISION(N)
    !                 the right hand side vector.
    ! 
    !         JOB     INTEGER
    !                 = 0         to solve  A*X = B ,
    !                 = nonzero   to solve  TRANS(A)*X = B , where
    !                             TRANS(A)  is the transpose.
    ! 
    !      On Return
    ! 
    !         B       the solution vector  X .
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
    ! ***ROUTINES CALLED  DAXPY, DDOT
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
    INTEGER LDA,N,ML,MU,IPVT(:),JOB
    DOUBLE PRECISION ABD(:,:),B(:)
    ! 
    DOUBLE PRECISION DDOT,T
    INTEGER K,KB,L,LA,LB,LM,M,NM1
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
       CALL DAXPY(LM,T,ABD(M+1:M+LM,K),1,B(K+1:K+LM),1)
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
       CALL DAXPY(LM,T,ABD(LA:LA + LM - 1,K),1,B(LB:LB + LM - 1),1)
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
       T = DDOT(LM,ABD(LA,K),1,B(LB),1)
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
       B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
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

#ifdef CUDA
  attributes(device) &
#endif  
  subroutine dgefa (a,lda,n,ipvt,info)

    !$acc routine seq
    !$acc routine(daxpy) seq
    !$acc routine(idamax) seq
    !$acc routine(dscal) seq

    integer lda,n,ipvt(:),info
    double precision a(:,:)
    ! 
    !      dgefa factors a double precision matrix by gaussian elimination.
    ! 
    !      dgefa is usually called by dgeco, but it can be called
    !      directly with a saving in time if  rcond  is not needed.
    !      (time for dgeco) = (1 + 9/n)*(time for dgefa) .
    ! 
    !      on entry
    ! 
    !         a       double precision(lda, n)
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
    !                 an integer vector of pivot indices.
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
    !      blas vdaxpy,dscal,idamax
    ! 
    !      internal variables
    ! 
    double precision t
    integer idamax,j,k,kp1,l,nm1
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
       l = idamax(n-k+1,a(k:n,k),1) + k - 1
       ipvt(k) = l
       ! 
       !         zero pivot implies this column already triangularized
       ! 
       if (a(l,k) .eq. 0.0d0) goto 40
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
       t = -1.0d0/a(k,k)
       call dscal(n-k,t,a(k+1:n,k),1)
       ! 
       !            row elimination with column indexing
       ! 
       do j = kp1, n
          t = a(l,j)
          if (l .eq. k) goto 20
          a(l,j) = a(k,j)
          a(k,j) = t
20        continue
          call daxpy(n-k,t,a(k+1:n,k),1,a(k+1:n,j),1)
       enddo
       goto 50
40     continue
       info = k
50     continue
    enddo
70  continue
    ipvt(n) = n
    if (a(n,n) .eq. 0.0d0) info = n
  end subroutine dgefa

#ifdef CUDA
  attributes(device) &
#endif  
  double precision function vddot (n,dx,incx,dy,incy)

    !$acc routine seq

    ! 
    !      forms the dot product of two vectors.
    !      uses unrolled loops for increments equal to one.
    !      jack dongarra, linpack, 3/11/78.
    ! 
    double precision dx(:),dy(:),dtemp
    integer i,incx,incy,ix,iy,m,mp1,n
    ! 
    vddot = 0.0d0
    dtemp = 0.0d0
    if (n.le.0) return
    if (incx.eq.1.and.incy.eq.1) go to 20
    ! 
    !      code for unequal increments or equal increments not equal to 1
    ! 
    ix = 1
    iy = 1
    if (incx.lt.0) ix = (-n+1)*incx + 1
    if (incy.lt.0) iy = (-n+1)*incy + 1
    do i = 1,n
       dtemp = dtemp + dx(ix)*dy(iy)
       ix = ix + incx
       iy = iy + incy
    enddo
    vddot = dtemp
    return
    ! 
    !      code for both increments equal to 1
    ! 
    ! 
    !      clean-up loop
    ! 
20  m = mod(n,5)
    if ( m .eq. 0 ) go to 40
    do i = 1,m
       dtemp = dtemp + dx(i)*dy(i)
    enddo
    if( n .lt. 5 ) go to 60
40  mp1 = m + 1
    do i = mp1,n,5
       dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
            dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
    enddo
60  vddot = dtemp
  end function vddot
  
end module linpack_module
