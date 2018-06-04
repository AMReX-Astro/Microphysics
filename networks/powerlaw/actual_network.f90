module actual_network

  use bl_types

  implicit none

  integer, parameter :: nspec = 3
  integer, parameter :: nspec_evolve = 2
  integer, parameter :: naux  = 0

  integer, parameter :: ifuel_  = 1
  integer, parameter :: iash_   = 2
  integer, parameter :: iinert_ = 3

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  double precision, save :: aion(nspec), zion(nspec), ebin(nspec)

  integer, parameter :: nrates = 0
  integer, parameter :: num_rate_groups = 0

contains
  
  subroutine actual_network_init

    use extern_probin_module, only: specific_q_burn
    use bl_constants_module, only: ZERO
    use fundamental_constants_module, only: N_A

    spec_names(ifuel_)  = "fuel"
    spec_names(iash_)   = "ash"
    spec_names(iinert_) = "inert"

    short_spec_names(ifuel_)  = "fuel"
    short_spec_names(iash_)   = "ash"
    short_spec_names(iinert_) = "inert"

    ! we are modeling a reaction f + f -> a + gamma, so baryon and
    ! charge conservation require that A_f = A_a / 2 and Z_f = Z_a / 2
    ! inert is a third species that doesn't participate in any 
    ! reactions
    aion(ifuel_)  = 2.0d0
    aion(iash_)   = 4.0d0
    aion(iinert_) = 8.0d0
    
    zion(ifuel_)  = 1.0d0
    zion(iash_)   = 2.0d0
    zion(iinert_) = 4.0d0

    ! Binding energies in erg / g

    ebin(ifuel_)  = ZERO
    ebin(iash_)   = specific_q_burn
    ebin(iinert_) = ZERO

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    ! Nothing to do here.

  end subroutine actual_network_finalize

end module actual_network
