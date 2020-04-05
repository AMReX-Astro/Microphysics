module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: ih1   = 1
  integer, parameter :: ihe3  = 2
  integer, parameter :: ihe4  = 3
  integer, parameter :: ic12  = 4
  integer, parameter :: in14  = 5
  integer, parameter :: io16  = 6
  integer, parameter :: ine20 = 7
  integer, parameter :: img24 = 8
  integer, parameter :: isi28 = 9
  integer, parameter :: is32  = 10
  integer, parameter :: iar36 = 11
  integer, parameter :: ica40 = 12
  integer, parameter :: iti44 = 13
  integer, parameter :: icr48 = 14
  integer, parameter :: icr56 = 15
  integer, parameter :: ife52 = 16
  integer, parameter :: ife54 = 17
  integer, parameter :: ife56 = 18
  integer, parameter :: ini56 = 19
  integer, parameter :: ineut = 20
  integer, parameter :: iprot = 21

  real(rt)        , allocatable :: bion(:), mion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion, mion
#endif

  character (len=32), parameter :: network_name = "aprox21"

  ! Some fundamental physical constants

  real(rt)        , parameter :: avo = 6.0221417930e23_rt
  real(rt)        , parameter :: c_light = 2.99792458e10_rt

  real(rt)        , parameter :: ev2erg  = 1.60217648740e-12_rt
  real(rt)        , parameter :: mev2erg = ev2erg*1.0e6_rt
  real(rt)        , parameter :: mev2gr  = mev2erg/c_light**2

  real(rt)        , parameter :: mn = 1.67492721184e-24_rt
  real(rt)        , parameter :: mp = 1.67262163783e-24_rt
  real(rt)        , parameter :: me = 9.1093821545e-28_rt

  ! Conversion factor for the nuclear energy generation rate.

  real(rt)        , parameter :: enuc_conv2 = -avo*c_light*c_light

  ! Rates data

  integer, parameter :: nrates = 112
  integer, parameter :: num_rate_groups = 4

  integer, parameter :: ir3a   = 1
  integer, parameter :: irg3a  = 2
  integer, parameter :: ircag  = 3
  integer, parameter :: ir1212 = 4
  integer, parameter :: ir1216 = 5
  integer, parameter :: ir1616 = 6
  integer, parameter :: iroga  = 7
  integer, parameter :: iroag  = 8
  integer, parameter :: irnega = 9
  integer, parameter :: irneag = 10
  integer, parameter :: irmgga = 11
  integer, parameter :: irmgag = 12
  integer, parameter :: irsiga = 13
  integer, parameter :: irmgap = 14
  integer, parameter :: iralpa = 15
  integer, parameter :: iralpg = 16
  integer, parameter :: irsigp = 17
  integer, parameter :: irsiag = 18
  integer, parameter :: irsga  = 19
  integer, parameter :: irsiap = 20
  integer, parameter :: irppa  = 21
  integer, parameter :: irppg  = 22
  integer, parameter :: irsgp  = 23
  integer, parameter :: irsag  = 24
  integer, parameter :: irarga = 25
  integer, parameter :: irsap  = 26
  integer, parameter :: irclpa = 27
  integer, parameter :: irclpg = 28
  integer, parameter :: irargp = 29
  integer, parameter :: irarag = 30
  integer, parameter :: ircaga = 31
  integer, parameter :: irarap = 32
  integer, parameter :: irkpa  = 33
  integer, parameter :: irkpg  = 34
  integer, parameter :: ircagp = 35
  integer, parameter :: ircaag = 36
  integer, parameter :: irtiga = 37
  integer, parameter :: ircaap = 38
  integer, parameter :: irscpa = 39
  integer, parameter :: irscpg = 40
  integer, parameter :: irtigp = 41
  integer, parameter :: irtiag = 42
  integer, parameter :: ircrga = 43
  integer, parameter :: irtiap = 44
  integer, parameter :: irvpa  = 45
  integer, parameter :: irvpg  = 46
  integer, parameter :: ircrgp = 47
  integer, parameter :: ircrag = 48
  integer, parameter :: irfega = 49
  integer, parameter :: ircrap = 50
  integer, parameter :: irmnpa = 51
  integer, parameter :: irmnpg = 52
  integer, parameter :: irfegp = 53
  integer, parameter :: irfeag = 54
  integer, parameter :: irniga = 55
  integer, parameter :: irfeap = 56
  integer, parameter :: ircopa = 57
  integer, parameter :: ircopg = 58
  integer, parameter :: irnigp = 59
  integer, parameter :: ir52ng = 60
  integer, parameter :: ir53gn = 61
  integer, parameter :: ir53ng = 62
  integer, parameter :: ir54gn = 63
  integer, parameter :: irfepg = 64
  integer, parameter :: ircogp = 65
  integer, parameter :: irheng = 66
  integer, parameter :: irhegn = 67
  integer, parameter :: irhng  = 68
  integer, parameter :: irdgn  = 69
  integer, parameter :: irdpg  = 70
  integer, parameter :: irhegp = 71
  integer, parameter :: irpen   = 72
  integer, parameter :: irnep   = 73
  integer, parameter :: irn56ec = 74
  integer, parameter :: irpp    = 75
  integer, parameter :: ir33    = 76
  integer, parameter :: irhe3ag = 77
  integer, parameter :: ircpg  = 78
  integer, parameter :: irnpg  = 79
  integer, parameter :: ifa    = 80
  integer, parameter :: ifg    = 81
  integer, parameter :: iropg  = 82
  integer, parameter :: irnag  = 83
  integer, parameter :: ir54ng   = 84
  integer, parameter :: ir55gn   = 85
  integer, parameter :: ir55ng   = 86
  integer, parameter :: ir56gn   = 87
  integer, parameter :: irfe54ap = 88
  integer, parameter :: irco57pa = 89
  integer, parameter :: irfe56pg = 90
  integer, parameter :: irco57gp = 91
  integer, parameter :: irr1   = 92
  integer, parameter :: irs1   = 93
  integer, parameter :: irt1   = 94
  integer, parameter :: iru1   = 95
  integer, parameter :: irv1   = 96
  integer, parameter :: irw1   = 97
  integer, parameter :: irx1   = 98
  integer, parameter :: ir1f54 = 99
  integer, parameter :: ir2f54 = 100
  integer, parameter :: ir3f54 = 101
  integer, parameter :: ir4f54 = 102
  integer, parameter :: ir5f54 = 103
  integer, parameter :: ir6f54 = 104
  integer, parameter :: ir7f54 = 105
  integer, parameter :: ir8f54 = 106
  integer, parameter :: iralf1 = 107
  integer, parameter :: iralf2 = 108
  integer, parameter :: irfe56_aux1 = 109
  integer, parameter :: irfe56_aux2 = 110
  integer, parameter :: irfe56_aux3 = 111
  integer, parameter :: irfe56_aux4 = 112

  character (len=20), save :: ratenames(nrates)

contains

  subroutine actual_network_init

    implicit none

    integer :: i

    call network_properties_init()

    allocate(bion(nspec))
    allocate(mion(nspec))
    

    ! Set the binding energy of the element
    bion(ih1)   = 0.0e0_rt
    bion(ihe3)  = 7.71819e0_rt
    bion(ihe4)  = 28.29603e0_rt
    bion(ic12)  = 92.16294e0_rt
    bion(in14)  = 104.65998e0_rt
    bion(io16)  = 127.62093e0_rt
    bion(ine20) = 160.64788e0_rt
    bion(img24) = 198.25790e0_rt
    bion(isi28) = 236.53790e0_rt
    bion(is32)  = 271.78250e0_rt
    bion(iar36) = 306.72020e0_rt
    bion(ica40) = 342.05680e0_rt
    bion(iti44) = 375.47720e0_rt
    bion(icr48) = 411.46900e0_rt
    bion(icr56) = 488.4970e0_rt
    bion(ife52) = 447.70800e0_rt
    bion(ife54) = 471.7696e0_rt
    bion(ife56) = 492.2450e0_rt
    bion(ini56) = 484.00300e0_rt
    bion(ineut) = 0.0e0_rt
    bion(iprot) = 0.0e0_rt

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

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

  end subroutine actual_network_init


  subroutine actual_network_finalize

    implicit none

    call network_properties_finalize()

    if (allocated(bion)) then
       deallocate(bion)
    endif
    if (allocated(mion)) then
       deallocate(mion)
    endif

  end subroutine actual_network_finalize

end module actual_network
