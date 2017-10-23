module dvode_dvnorm_module

  use dvode_constants_module

  implicit none

contains

#ifdef CUDA
  attributes(device) &
#endif  
  function dvnorm(N, V, W) result(dvn)

    !$acc routine seq
    
    ! ***BEGIN PROLOGUE  DVNORM
    ! ***SUBSIDIARY
    ! ***PURPOSE  Weighted root-mean-square vector norm.
    ! ***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This function routine computes the weighted root-mean-square norm
    !   of the vector of length N contained in the array V, with weights
    !   contained in the array W of length N:
    !     DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
    ! 
    ! ***SEE ALSO  DLSODE
    ! ***ROUTINES CALLED  (NONE)
    ! ***REVISION HISTORY  (YYMMDD)
    !    791129  DATE WRITTEN
    !    890501  Modified prologue to SLATEC/LDOC format.  (FNF)
    !    890503  Minor cosmetic changes.  (FNF)
    !    930809  Renamed to allow single/double precision versions. (ACH)
    ! ***END PROLOGUE  DVNORM
    ! **End

    use bl_types, only: dp_t

    implicit none

    integer    :: N, I
    real(dp_t) :: V(N), W(N)
    real(dp_t) :: SUM, dvn, dscratch

    SUM = 0.0D0
    do I = 1,N
       SUM = SUM + (V(I)*W(I))**2
    end do
    dscratch = SUM/N
    dvn = SQRT(dscratch)
    return
  end function dvnorm

end module dvode_dvnorm_module
