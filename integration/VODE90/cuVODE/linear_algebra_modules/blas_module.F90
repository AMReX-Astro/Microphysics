module blas_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
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
    !      forms the dot product of two arrays.
    !      uses unrolled loops for increments equal to one.
    !      jack dongarra, linpack, 3/11/78.
    !      modified 12/3/93, array(1) declarations changed to array(*)
    ! 
    ! 
    !      .. Local Scalars ..
    DOUBLE PRECISION DTEMP
    INTEGER I,IX,IY,M,MP1

    !$gpu

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

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  SUBROUTINE DSCALN(N,DA,DX,INCX)
    ! Only operates on arrays of size N

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
    !      scales a array by a constant.
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
    !$gpu
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

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  SUBROUTINE DSCAL(N,DA,DX,INCX)
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
    !      scales a array by a constant.
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
    !$gpu
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

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  FUNCTION IDAMAX(N,DX,INCX) result(index)
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
    INTEGER I,IX,index
    !      ..
    !      .. Intrinsic Functions ..
    INTRINSIC DABS
    !$gpu
    !      ..
    index = 0
    IF (N.LT.1 .OR. INCX.LE.0) RETURN
    index = 1
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
       index = I
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
       index = I
       DMAX = DABS(DX(I))
    end do
    RETURN
  END FUNCTION IDAMAX
  
end module blas_module
