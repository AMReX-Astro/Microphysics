module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

  implicit none

  integer, parameter :: nrates = 112

  character (len=20) :: ratenames(nrates)

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter, private :: avo = 6.0221417930d23
  double precision, parameter, private :: c_light = 2.99792458d10
  double precision, parameter, private :: enuc_conv2 = -avo*c_light*c_light

contains

  subroutine actual_burner(state_in, state_out, dt, time)

    use vode_module, only: vode_burner

    implicit none

    type (eos_t),        intent(in   ) :: state_in
    type (eos_t),        intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    call vode_burner(state_in, state_out, dt, time)

  end subroutine actual_burner



  subroutine actual_burner_init()

    use rates_module, only: rates_init
    use screening_module, only: screening_init
    use network_indices
    use rpar_indices

    implicit none

    ! set the names of the reaction rates
    ratenames(ir3a)   = 'r_he4_he4_he4_to_c12'
    ratenames(irg3a)  = 'r_c12_to_he4_he4_he4'
    ratenames(ircag)  = 'r_c12_ag_o16'
    ratenames(ir1212) = 'r1212'
    ratenames(ir1216) = 'r1216'
    ratenames(ir1616) = 'r1616'
    ratenames(iroga)  = 'r_o16_ga_c12'
    ratenames(iroag)  = 'r_o16_ag_ne20'
    ratenames(irnega) = 'r_ne20_ga_o16'
    ratenames(irneag) = 'r_ne20_ag_mg24'
    ratenames(irmgga) = 'r_mg24_ga_ne20'
    ratenames(irmgag) = 'r_mg24_ag_si28'
    ratenames(irsiga) = 'r_si28_ga_mg24'
    ratenames(irmgap) = 'r_mg24_ap_al27'
    ratenames(iralpa) = 'r_al27_pa_mg24'
    ratenames(iralpg) = 'r_al27_pg_si28'
    ratenames(irsigp) = 'r_si28_gp_al27'
    ratenames(irsiag) = 'r_si28_ag_s32'
    ratenames(irsga)  = 'r_s32_ga_si28'
    ratenames(irsiap) = 'r_si28_ap_p31'
    ratenames(irppa)  = 'r_p31_pa_si28'
    ratenames(irppg)  = 'r_p31_pg_s32'
    ratenames(irsgp)  = 'r_s32_gp_p31'
    ratenames(irsag)  = 'r_s32_ag_ar36'
    ratenames(irarga) = 'r_ar36_ga_s32'
    ratenames(irsap)  = 'r_s32_ap_cl35'
    ratenames(irclpa) = 'r_cl35_pa_s32'
    ratenames(irclpg) = 'r_cl35_pg_ar36'
    ratenames(irargp) = 'r_ar36_gp_cl35'
    ratenames(irarag) = 'r_ar36_ag_ca40'
    ratenames(ircaga) = 'r_ca40_ga_ar36'
    ratenames(irarap) = 'r_ar36_ap_k39'
    ratenames(irkpa)  = 'r_k39_pa_ar36'
    ratenames(irkpg)  = 'r_k39_pg_ca40'
    ratenames(ircagp) = 'r_ca40_gp_k39'
    ratenames(ircaag) = 'r_ca40_ag_ti44'
    ratenames(irtiga) = 'r_ti44_ga_ca40'
    ratenames(ircaap) = 'r_ca40_ap_sc43'
    ratenames(irscpa) = 'r_sc43_pa_ca40'
    ratenames(irscpg) = 'r_sc43_pg_ti44'
    ratenames(irtigp) = 'r_ti44_gp_sc43'
    ratenames(irtiag) = 'r_ti44_ag_cr48'
    ratenames(ircrga) = 'r_cr48_ga_ti44'
    ratenames(irtiap) = 'r_ti44_ap_v47'
    ratenames(irvpa)  = 'r_v47_pa_ti44'
    ratenames(irvpg)  = 'r_v47_pg_cr48'
    ratenames(ircrgp) = 'r_cr48_gp_v47'
    ratenames(ircrag) = 'r_cr48_ag_fe52'
    ratenames(irfega) = 'r_fe52_ga_cr48'
    ratenames(ircrap) = 'r_cr48_ap_mn51'
    ratenames(irmnpa) = 'r_mn51_pa_cr48'
    ratenames(irmnpg) = 'r_mn51_pg_fe52'
    ratenames(irfegp) = 'r_fe52_gp_mn51'
    ratenames(irfeag) = 'r_fe52_ag_ni56'
    ratenames(irniga) = 'r_ni56_ga_fe52'
    ratenames(irfeap) = 'r_fe52_ap_co55'
    ratenames(ircopa) = 'r_co55_pa_fe52'
    ratenames(ircopg) = 'r_co55_pg_ni56'
    ratenames(irnigp) = 'r_ni56_gp_co55'

    ! for fe54 photodisintegration
    ratenames(ir52ng) = 'r_fe52_ng_fe53'
    ratenames(ir53gn) = 'r_fe53_gn_fe52'
    ratenames(ir53ng) = 'r_fe53_ng_fe54'
    ratenames(ir54gn) = 'r_fe54_gn_fe53'
    ratenames(irfepg) = 'r_fe54_pg_co55'
    ratenames(ircogp) = 'r_co55_gp_fe54'

    ! for he4 photodisintegration
    ratenames(irheng)  = 'r_he3_ng_he4'
    ratenames(irhegn)  = 'r_he4_gn_he3'
    ratenames(irhng)   = 'r_h1_ng_h2'
    ratenames(irdgn)   = 'r_h2_gn_h1'
    ratenames(irdpg)   = 'r_h2_pg_he3'
    ratenames(irhegp)  = 'r_he3_gp_h2'

    ! for weak reactions
    ratenames(irpen)   = 'r_prot_to_neut'
    ratenames(irnep)   = 'r_neut_to_prot'
    ratenames(irn56ec) = 'r_ni56ec_to_fe56'

    ! ppchain
    ratenames(irpp)    = 'rpp_to_he3'
    ratenames(ir33)    = 'r_he3_he3_to_h1_h1_he4'
    ratenames(irhe3ag) = 'r_he3_ag_be7'

    ! cno cycles
    ratenames(ircpg)   = 'r_c12_pg_n13'
    ratenames(irnpg)   = 'r_n14_pg_o15'
    ratenames(iropg)   = 'r_o16_pg_f17'
    ratenames(ifa)     = 'rfa'
    ratenames(ifg)     = 'rfg'
    ratenames(irnag)   = 'r_n14_ag_f18'

    ! for reactions to fe56 
    ratenames(ir54ng)   = 'r_fe54_ng_fe55'
    ratenames(ir55gn)   = 'r_fe55_gn_fe54'
    ratenames(ir55ng)   = 'r_fe55_ng_fe56'
    ratenames(ir56gn)   = 'r_fe56_gn_fe55'
    ratenames(irfe54ap) = 'r_fe54_ap_co57'
    ratenames(irco57pa) = 'r_co57_pa_fe54'
    ratenames(irfe56pg) = 'r_fe56_pg_co57'
    ratenames(irco57gp) = 'r_co57_gp_fe56'


    ! the equilibrium links
    ratenames(irr1)   = 'r_al27_equil'
    ratenames(irs1)   = 'r_p31_equil'
    ratenames(irt1)   = 'r_cl35_equil'
    ratenames(iru1)   = 'r_k39_equil'
    ratenames(irv1)   = 'r_sc43_equil'
    ratenames(irw1)   = 'r_v47_equil'
    ratenames(irx1)   = 'r_mn51_equil'

    ratenames(ir1f54) = 'r1f54'
    ratenames(ir2f54) = 'r2f54'
    ratenames(ir3f54) = 'r3f54'
    ratenames(ir4f54) = 'r4f54'
    ratenames(ir5f54) = 'r5f54'
    ratenames(ir6f54) = 'r6f54'
    ratenames(ir7f54) = 'r7f54'
    ratenames(ir8f54) = 'r8f54'

    ratenames(iralf1) = 'ralf1'
    ratenames(iralf2) = 'ralf2'

    ratenames(irfe56_aux1) = 'rfe56aux1'
    ratenames(irfe56_aux2) = 'rfe56aux2'
    ratenames(irfe56_aux3) = 'rfe56aux3'
    ratenames(irfe56_aux4) = 'rfe56aux4'

    call init_rpar_indices(nrates, nspec)

    ! Add data to store the rates

    irp_dydt  = get_next_rpar_index(nspec)
    irp_rates = get_next_rpar_index(nrates)
    irp_drdy1 = get_next_rpar_index(nrates)
    irp_drdy2 = get_next_rpar_index(nrates)
    
    call rates_init()

    call set_up_screening_factors()

    call screening_init()

  end subroutine actual_burner_init



  ! Compute and store the more expensive screening factors  

  subroutine set_up_screening_factors()

    use screening_module, only: add_screening_factor
    use network, only: aion, zion
    use network_indices

    implicit none

    call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ihe4),aion(ihe4),4.0d0,8.0d0)
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))
    
    call add_screening_factor(zion(io16),aion(io16),zion(io16),aion(io16))
    
    call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(13.0d0,27.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(isi28),aion(isi28),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(15.0d0,31.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(is32),aion(is32),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(17.0d0,35.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(iar36),aion(iar36),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(19.0d0,39.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(ica40),aion(ica40),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(21.0d0,43.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(iti44),aion(iti44),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(23.0d0,47.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(icr48),aion(icr48),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(25.0d0,51.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(ife52),aion(ife52),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(27.0d0,55.0d0,1.0d0,1.0d0)
    
    call add_screening_factor(zion(ife54),aion(ife54),1.0d0,1.0d0)
    
    call add_screening_factor(zion(ife54),aion(ife54),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ife56),aion(ife56),1.0d0,1.0d0)
    
    call add_screening_factor(1.0d0,2.0d0,zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(ih1),aion(ih1),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3))
    
    call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4))
    
    call add_screening_factor(zion(ic12),aion(ic12),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(in14),aion(in14),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(io16),aion(io16),zion(ih1),aion(ih1))
    
    call add_screening_factor(zion(in14),aion(in14),zion(ihe4),aion(ihe4))

  end subroutine set_up_screening_factors



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    implicit none

    double precision :: dydt(nspec), enuc

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_burner_module
