module cuvode_dvnorm_module

  use amrex_fort_module, only: rt => amrex_real
  use cuvode_parameters_module, only: VODE_LMAX, VODE_NEQS, VODE_LIW,   &
                                      VODE_LENWM, VODE_MAXORD

  use cuvode_constants_module

  implicit none

contains

#if defined(AMREX_USE_CUDA) && !defined(AMREX_USE_GPU_PRAGMA)
  attributes(device) &
#endif
  function dvnorm(V, W) result(dvn)

    !$acc routine seq
    
    ! ***BEGIN PROLOGUE  DVNORM
    ! ***SUBSIDIARY
    ! ***PURPOSE  Weighted root-mean-square array norm.
    ! ***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
    ! ***AUTHOR  Hindmarsh, Alan C., (LLNL)
    ! ***DESCRIPTION
    ! 
    !   This function routine computes the weighted root-mean-square norm
    !   of the array of length N contained in the array V, with weights
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

    implicit none

    ! Declare arguments
    real(rt), intent(in   ) :: V(VODE_NEQS), W(VODE_NEQS)

    ! Declare return variable
    real(rt) :: dvn

    ! Declare local variables
    real(rt) :: SUM
    integer    :: I

    !$gpu

    SUM = 0.0e0_rt
    do I = 1,VODE_NEQS
       SUM = SUM + (V(I)*W(I))**2
    end do
    dvn = SQRT(SUM/VODE_NEQS)
    RETURN
  end function dvnorm

end module cuvode_dvnorm_module
