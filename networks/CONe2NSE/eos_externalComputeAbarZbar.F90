!
! Dean Townsley 2008
!
!  This function is a callback for the Eos/Helmholtz/ExternelAbarZbar
!  EOS implementation.  (e.g. the interface is declared there)
!  the local Abar and Zbar are calculated from the provided mass scalars

subroutine eos_externalComputeAbarZbar(solnScalars, abarData, zbarData)

  use bn_paraInterface, ONLY: bn_paraFuelAshProperties

  implicit none

#include "constants.h"
#include "Flash.h"

  real, intent(in), dimension(:,:)  :: solnScalars
  real, intent(out), dimension(:)                      :: abarData, zbarData

  integer :: c12ii = CI_MSCALAR-SPECIES_BEGIN+1
  integer :: ne22ii = NEI_MSCALAR-SPECIES_BEGIN+1
  integer :: phifai = PHFA_MSCALAR-SPECIES_BEGIN+1
  integer :: phiaqi = PHAQ_MSCALAR-SPECIES_BEGIN+1
  integer :: phiqni = PHQN_MSCALAR-SPECIES_BEGIN+1
  integer :: yei = YE_MSCALAR-SPECIES_BEGIN+1
  integer :: dyiqni = DYQN_MSCALAR-SPECIES_BEGIN+1

  real :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a
  real :: phi_fa, phi_aq, dyi_qn
  integer :: i

  do i = 1, ubound(solnScalars,2)
     call bn_paraFuelAshProperties( solnScalars(c12ii,i), &
                                    solnScalars(ne22ii,i), &
                                    ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a )
     ! truncate advection errors just as is done in bn_paraBurn
     phi_fa = max(0.0,min(1.0,solnScalars(phifai,i)))
     phi_aq = max(0.0,min(phi_fa,solnScalars(phiaqi,i)))
     dyi_qn = max(0.0,min(1.0,solnScalars(dyiqni,i)))
     abarData(i) = 1.0 / ( (1.0-phi_fa)*yi_f + (phi_fa-phi_aq)*yi_a + dyi_qn  )
     zbarData(i) = abarData(i) * solnScalars(yei,i)

     ! check
     if ( .not. (abarData(i) > 0.0 ) ) then
        !print but don't die, EOS will die and report density and temperature
        write (6,*) "[eos_externalComputeAbarZbar] negative abar =", abarData(i)
        write (6,*) "grid data ci, nei, phfa, phaq, phqn, ye, dyqn :"
        write (6,*) solnScalars(c12ii,i), solnScalars(ne22ii,i), solnScalars(phifai,i), &
                    solnScalars(phiaqi,i), solnScalars(phiqni,i), solnScalars(yei,i), &
                    solnScalars(dyiqni,i)
     endif
  enddo

end subroutine
