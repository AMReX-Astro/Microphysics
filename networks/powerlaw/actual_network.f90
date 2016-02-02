module actual_network

  use bl_types
  use actual_network_data

  implicit none

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

end module actual_network
