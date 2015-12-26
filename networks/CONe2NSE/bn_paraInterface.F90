! interfaces for private functions for parametric Burn unit
!
! Dean Townsley 2008
!

module bn_paraInterface

#include "constants.h"
#include "Flash.h" 

  interface bn_paraBurn
     subroutine bn_paraBurn(dens, temp, eint, pres, &
                       xc12init, xne22init, &
                       flame, flamedot, phi_fa, phi_aq, phi_qn, &
                       ye, dyi_qn, dqbar_qn, &
                       qdot, edotnu, dt, react_proximity, shock, &
                       ignite_detonation, phi_fa_det)
        implicit none
        real, intent(in)    :: dens, temp, eint, pres
        real, intent(in)    :: xc12init, xne22init
        real, intent(in)    :: flame, flamedot
        real, intent(inout) :: phi_fa, phi_aq, phi_qn
        real, intent(inout) :: ye, dyi_qn, dqbar_qn
        real, intent(out)   :: qdot, edotnu
        real, intent(in)    :: dt, react_proximity, shock
        logical, intent(in) :: ignite_detonation
        real, intent(in)    :: phi_fa_det
     end subroutine bn_paraBurn
  end interface
  
  interface bn_paraFuelAshProperties
     subroutine bn_paraFuelAshProperties(xc12initial, xne22initial, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)
        implicit none
        real, intent(in)  :: xc12initial, xne22initial
        real, intent(out) :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a
     end subroutine bn_paraFuelAshProperties
  end interface

  interface bn_paraAllIgnite
     subroutine bn_paraAllIgnite(ignition_conditions, det_num)
        implicit none
        integer, intent(in)  :: det_num
        logical, allocatable, dimension(:), intent(inout) :: ignition_conditions
     end subroutine bn_paraAllIgnite
  end interface

  interface bn_paraSpark
     subroutine bn_paraSpark(x, y, z, dens, pres, phfa, flame, c12, ne22, &
                             detX, detY, detZ, ignition_test)
        implicit none
        real, intent(in) :: x, y, z, dens, pres, phfa, flame, c12, ne22
        real, intent(out) :: detX, detY, detZ
        logical, intent(out) :: ignition_test
     end subroutine bn_paraSpark
  end interface

  interface bn_paraAllSpark
     subroutine bn_paraAllSpark(ignition_coords, det_num, &
                                det_xCoord, det_yCoord, det_zCoord)
        use Burn_data, ONLY : pbIgnNumMax
        implicit none
        real, allocatable, dimension(:), intent(inout) :: ignition_coords
        integer, intent(inout) :: det_num
        real, allocatable, dimension(:), intent(out) :: det_xCoord, &
                                                        det_yCoord, &
                                                        det_zCoord
     end subroutine bn_paraAllSpark
  end interface

end module bn_paraInterface
