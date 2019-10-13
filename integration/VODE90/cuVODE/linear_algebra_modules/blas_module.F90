module blas_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  function idamax(N,DX) result(index)
    !      .. Scalar Arguments ..
    INTEGER N
    !      ..
    !      .. Array Arguments ..
    DOUBLE PRECISION DX(N)
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

    !$gpu

    index = 1
    DMAX = DABS(DX(1))
    DO I = 2,N
       IF (DABS(DX(I)).LE.DMAX) cycle
       index = I
       DMAX = DABS(DX(I))
    end do
  end function idamax

end module blas_module
