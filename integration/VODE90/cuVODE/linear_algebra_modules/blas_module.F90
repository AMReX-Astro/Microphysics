module blas_module

  implicit none

contains

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
