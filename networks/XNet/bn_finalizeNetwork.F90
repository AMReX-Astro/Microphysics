!!****f* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/bn_finalizeNetwork
!!
!! NAME
!!  
!!  bn_finalizeNetwork
!!
!!
!! SYNOPSIS
!! 
!!  call bn_finalizeNetwork()
!!
!!  
!! DESCRIPTION
!!
!!  Finalizes the bn_* submodule
!!
!!***


subroutine bn_finalizeNetwork()

  use bn_interface, ONLY : bn_xnetFinalize

  implicit none

  call bn_xnetFinalize

  return

end subroutine bn_finalizeNetwork
