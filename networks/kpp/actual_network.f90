module actual_network

  use bl_types

  implicit none

  integer, parameter :: nspec = 2
  integer, parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)

  integer, parameter :: ifuel_ = 1
  integer, parameter :: iash_ = 2

contains
  
  subroutine actual_network_init

    spec_names(ifuel_)  = "fuel"
    spec_names(iash_)  = "ash"

    short_spec_names(ifuel_)  = "fuel"
    short_spec_names(iash_)  = "ash"

    aion(ifuel_)  = 2.0_dp_t
    aion(iash_)  = 4.0_dp_t
    
    zion(ifuel_)  = 1.0_dp_t
    zion(iash_)  = 2.0_dp_t

  end subroutine actual_network_init

end module actual_network
