module actual_network

  use amrex_fort_module, only : rt => amrex_real

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: nspec = 2
  integer, parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), save :: aion(nspec), zion(nspec), ebin(nspec)

  integer, parameter :: ifuel_ = 1
  integer, parameter :: iash_ = 2

contains
  
  subroutine actual_network_init

    spec_names(ifuel_)  = "fuel"
    spec_names(iash_)  = "ash"

    short_spec_names(ifuel_)  = "fuel"
    short_spec_names(iash_)  = "ash"

    aion(ifuel_)  = 2.0_rt
    aion(iash_)  = 4.0_rt
    
    zion(ifuel_)  = 1.0_rt
    zion(iash_)  = 2.0_rt

  end subroutine actual_network_init

end module actual_network
end module actual_network
