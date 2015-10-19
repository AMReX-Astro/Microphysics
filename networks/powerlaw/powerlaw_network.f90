module actual_network

  use bl_types

  implicit none

  integer, parameter :: nspec = 2
  integer, parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec) 
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save :: aion(nspec), zion(nspec)

  integer, parameter :: ifuel_ = 1
  integer, parameter :: iash_ = 2

contains
  
  subroutine actual_network_init

    use rpar_indices

    spec_names(ifuel_)  = "fuel"
    spec_names(iash_)  = "ash"

    short_spec_names(ifuel_)  = "fuel"
    short_spec_names(iash_)  = "ash"

    ! we are modeling a reaction f + f -> a + gamma, so baryon and
    ! charge conservation require that A_f = A_a / 2 and Z_f = Z_a / 2
    aion(ifuel_)  = 2.0d0
    aion(iash_)  = 4.0d0
    
    zion(ifuel_)  = 1.0d0
    zion(iash_)  = 2.0d0

    call init_rpar_indices(nspec)

  end subroutine actual_network_init

end module actual_network
