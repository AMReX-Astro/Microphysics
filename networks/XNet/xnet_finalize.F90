subroutine xnet_finalize()

   use xnet_interface, ONLY : jacobian_finalize

   implicit none

   call jacobian_finalize

   return
end subroutine bn_xnetFinalize
