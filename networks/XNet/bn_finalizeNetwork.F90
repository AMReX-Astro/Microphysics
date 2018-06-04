subroutine bn_finalizeNetwork()

  use bn_interface, ONLY : bn_xnetFinalize

  implicit none

  call bn_xnetFinalize

  return

end subroutine bn_finalizeNetwork
