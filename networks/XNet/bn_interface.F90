!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/bn_interface
!!
!! SYNOPSIS
!!   use bn_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Burn module that defines its
!! private interfaces.
!!
!!***

Module bn_interface

#include "Flash.h"
#include "constants.h"

  interface
     subroutine bn_xnetInit(data_dir,data_desc)
        implicit none
        character(*), intent(in) :: data_dir
        character(80), intent(out) :: data_desc
     end subroutine bn_xnetInit
  end interface

  interface
     subroutine bn_xnetFinalize()
        implicit none
     end subroutine bn_xnetFinalize
  end interface

  interface
     subroutine bn_burner(tstep,temp,density,xIn,xOut,sdotRate,burnedZone,kstep)
       implicit none
       logical, intent(IN), dimension(:)              :: burnedZone
       real, intent(IN)                               :: tstep
       real, intent(IN), dimension(:)                 :: temp,density
       real, intent(OUT), dimension(size(temp))       :: sdotRate
       real, intent(IN), dimension(:,:)               :: xIn
       real, intent(OUT), dimension(size(xIn,1),size(xIn,2)) :: xOut
       integer, intent(OUT)                           :: kstep
     end subroutine bn_burner
  end interface

  interface
     subroutine bn_azbar()
       implicit none
     end subroutine bn_azbar
  end interface

  interface
     subroutine bn_ecapnuc(etakep,temp,rpen,rnep,spen,snep)
       implicit none
       real, intent(IN)     :: etakep, temp
       real, intent(OUT)    :: rpen,rnep,spen,snep
     end subroutine bn_ecapnuc
  end interface

  interface
     real function bn_ifermi12(f)
       implicit none
       real, intent(IN) :: f
     end function bn_ifermi12
  end interface

  interface
     subroutine bn_mazurek(btemp,bden,y56,ye,rn56ec,sn56ec)
       implicit none
       real, intent(IN) :: btemp, bden
       real, intent(IN) :: y56, ye
       real, intent(OUT):: rn56ec,sn56ec

     end subroutine bn_mazurek
  end interface

  interface
     subroutine bn_mcord(i,j,iloc,jloc,nzo,np,eloc,nterm,np2)
       implicit none
       integer, intent(IN)  ::  i,j, np, np2
       integer, intent(INOUT) :: nterm, nzo
       integer, intent(OUT) ::  iloc(np),jloc(np),eloc(np2)
     end subroutine bn_mcord
  end interface

  interface
     subroutine bn_screen4(zbarr,abarr,z2barr,z1,a1,z2,a2, & 
          &                   jscreen,init,scorr,scorrdt)
       implicit none
       integer, intent(IN)   :: jscreen, init
       real, intent(IN)      :: abarr, zbarr, z2barr, z1, a1, z2, a2
       real, intent(OUT)     :: scorr
       real, intent(OUT), optional     :: scorrdt
     end subroutine bn_screen4
  end interface

  interface
     subroutine bn_sneutx()
       implicit none
     end subroutine bn_sneutx
  end interface

end Module bn_interface
