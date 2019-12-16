module actual_network

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, parameter :: nspec = 19
  integer, parameter :: nspec_evolve = 19
  integer, parameter :: naux  = 0

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
  integer, parameter :: ife52 = 15
  integer, parameter :: ife54 = 16
  integer, parameter :: ini56 = 17
  integer, parameter :: ineut = 18
  integer, parameter :: iprot = 19

  real(rt)        , allocatable :: aion(:), zion(:), nion(:)
  real(rt)        , allocatable :: bion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, nion, bion, mion, wion
#endif

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  character (len=32), parameter :: network_name = "aprox19"

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

  integer, parameter :: nrates = 100
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
  integer, parameter :: irr1   = 84
  integer, parameter :: irs1   = 85
  integer, parameter :: irt1   = 86
  integer, parameter :: iru1   = 87
  integer, parameter :: irv1   = 88
  integer, parameter :: irw1   = 89
  integer, parameter :: irx1   = 90
  integer, parameter :: ir1f54 = 91
  integer, parameter :: ir2f54 = 92
  integer, parameter :: ir3f54 = 93
  integer, parameter :: ir4f54 = 94
  integer, parameter :: ir5f54 = 95
  integer, parameter :: ir6f54 = 96
  integer, parameter :: ir7f54 = 97
  integer, parameter :: ir8f54 = 98

  integer, parameter :: iralf1 = 99
  integer, parameter :: iralf2 = 100

  character (len=16), save :: ratenames(nrates)
  
contains
  
  subroutine actual_network_init

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    ! The following comes directly from init_aprox19

    short_spec_names(ih1)   = 'h1'
    short_spec_names(ihe3)  = 'he3'  
    short_spec_names(ihe4)  = 'he4'
    short_spec_names(ic12)  = 'c12'
    short_spec_names(in14)  = 'n14'
    short_spec_names(io16)  = 'o16'
    short_spec_names(ine20) = 'ne20'
    short_spec_names(img24) = 'mg24'
    short_spec_names(isi28) = 'si28'
    short_spec_names(is32)  = 's32'
    short_spec_names(iar36) = 'ar36'
    short_spec_names(ica40) = 'ca40'
    short_spec_names(iti44) = 'ti44'
    short_spec_names(icr48) = 'cr48'
    short_spec_names(ife52) = 'fe52'
    short_spec_names(ife54) = 'fe54'
    short_spec_names(ini56) = 'ni56'
    short_spec_names(ineut) = 'neut'
    short_spec_names(iprot) = 'prot'

    spec_names(ih1)   = "hydrogen-1"
    spec_names(ihe3)  = "helium-3"
    spec_names(ihe4)  = "helium-4"
    spec_names(ic12)  = "carbon-12"
    spec_names(in14)  = "nitrogen-14"
    spec_names(io16)  = "oxygen-16"
    spec_names(ine20) = "neon-20"
    spec_names(img24) = "magnesium-24"
    spec_names(isi28) = "silicon-28"
    spec_names(is32)  = "sulfur-32"
    spec_names(iar36) = "argon-36"
    spec_names(ica40) = "calcium-40"
    spec_names(iti44) = "titanium-44"
    spec_names(icr48) = "chromium-48"
    spec_names(ife52) = "iron-52"
    spec_names(ife54) = "iron-54"
    spec_names(ini56) = "nickel-56"
    spec_names(ineut) = "neutron"
    spec_names(iprot) = "proton"
    
    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(nion(nspec))
    allocate(bion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    ! Set the number of nucleons in the element
    aion(ih1)   = 1.0e0_rt
    aion(ihe3)  = 3.0e0_rt
    aion(ihe4)  = 4.0e0_rt
    aion(ic12)  = 12.0e0_rt
    aion(in14)  = 14.0e0_rt
    aion(io16)  = 16.0e0_rt
    aion(ine20) = 20.0e0_rt
    aion(img24) = 24.0e0_rt
    aion(isi28) = 28.0e0_rt
    aion(is32)  = 32.0e0_rt
    aion(iar36) = 36.0e0_rt
    aion(ica40) = 40.0e0_rt
    aion(iti44) = 44.0e0_rt
    aion(icr48) = 48.0e0_rt
    aion(ife52) = 52.0e0_rt
    aion(ife54) = 54.0e0_rt
    aion(ini56) = 56.0e0_rt
    aion(ineut) = 1.0e0_rt
    aion(iprot) = 1.0e0_rt

    ! Set the number of protons in the element
    zion(ih1)   = 1.0e0_rt
    zion(ihe3)  = 2.0e0_rt
    zion(ihe4)  = 2.0e0_rt
    zion(ic12)  = 6.0e0_rt
    zion(in14)  = 7.0e0_rt
    zion(io16)  = 8.0e0_rt
    zion(ine20) = 10.0e0_rt
    zion(img24) = 12.0e0_rt
    zion(isi28) = 14.0e0_rt
    zion(is32)  = 16.0e0_rt
    zion(iar36) = 18.0e0_rt
    zion(ica40) = 20.0e0_rt
    zion(iti44) = 22.0e0_rt
    zion(icr48) = 24.0e0_rt
    zion(ife52) = 26.0e0_rt
    zion(ife54) = 26.0e0_rt
    zion(ini56) = 28.0e0_rt
    zion(ineut) = 0.0e0_rt
    zion(iprot) = 1.0e0_rt    

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
    bion(ife52) = 447.70800e0_rt
    bion(ife54) = 471.7696e0_rt
    bion(ini56) = 484.00300e0_rt
    bion(ineut) = 0.0e0_rt
    bion(iprot) = 0.0e0_rt

    ! Set the number of neutrons
    nion(:) = aion(:) - zion(:)

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    wion(:) = aion(:)

    ! set the names of the reaction rates
    ratenames(ir3a)   = 'r3a  '
    ratenames(irg3a)  = 'rg3a '
    ratenames(ircag)  = 'rcag '
    ratenames(ir1212) = 'r1212'
    ratenames(ir1216) = 'r1216'
    ratenames(ir1616) = 'r1616'
    ratenames(iroga)  = 'roga '
    ratenames(iroag)  = 'roag '
    ratenames(irnega) = 'rnega'
    ratenames(irneag) = 'rneag'
    ratenames(irmgga) = 'rmgga'
    ratenames(irmgag) = 'rmgag'
    ratenames(irsiga) = 'rsiga'
    ratenames(irmgap) = 'rmgap'
    ratenames(iralpa) = 'ralpa'
    ratenames(iralpg) = 'ralpg'
    ratenames(irsigp) = 'rsigp'
    ratenames(irsiag) = 'rsiag'
    ratenames(irsga)  = 'rsga '
    ratenames(irsiap) = 'rsiap'
    ratenames(irppa)  = 'rppa '
    ratenames(irppg)  = 'rppg '
    ratenames(irsgp)  = 'rsgp '
    ratenames(irsag)  = 'rsag '
    ratenames(irarga) = 'rarga'
    ratenames(irsap)  = 'rsap '
    ratenames(irclpa) = 'rclpa'
    ratenames(irclpg) = 'rclpg'
    ratenames(irargp) = 'rargp'
    ratenames(irarag) = 'rarag'
    ratenames(ircaga) = 'rcaga'
    ratenames(irarap) = 'rarap'
    ratenames(irkpa)  = 'rkpa '
    ratenames(irkpg)  = 'rkpg '
    ratenames(ircagp) = 'rcagp'
    ratenames(ircaag) = 'rcaag'
    ratenames(irtiga) = 'rtiga'
    ratenames(ircaap) = 'rcaap'
    ratenames(irscpa) = 'rscpa'
    ratenames(irscpg) = 'rscpg'
    ratenames(irtigp) = 'rtigp'
    ratenames(irtiag) = 'rtiag'
    ratenames(ircrga) = 'rcrga'
    ratenames(irtiap) = 'rtiap'
    ratenames(irvpa)  = 'rvpa '
    ratenames(irvpg)  = 'rvpg '
    ratenames(ircrgp) = 'rcrgp'
    ratenames(ircrag) = 'rcrag'
    ratenames(irfega) = 'rfega'
    ratenames(ircrap) = 'rcrap'
    ratenames(irmnpa) = 'rmnpa'
    ratenames(irmnpg) = 'rmnpg'
    ratenames(irfegp) = 'rfegp'
    ratenames(irfeag) = 'rfeag'
    ratenames(irniga) = 'rniga'
    ratenames(irfeap) = 'rfeap'
    ratenames(ircopa) = 'rcopa'
    ratenames(ircopg) = 'rcopg'
    ratenames(irnigp) = 'rnigp'

    ! for fe54 photodisintegration
    ratenames(ir52ng) = 'r52ng'
    ratenames(ir53gn) = 'r53gn'
    ratenames(ir53ng) = 'r53ng'
    ratenames(ir54gn) = 'r54gn'
    ratenames(irfepg) = 'rfepg'
    ratenames(ircogp) = 'rcogp'

    ! for he4 photodisintegration
    ratenames(irheng)  = 'rheng'
    ratenames(irhegn)  = 'rhegn'
    ratenames(irhng)   = 'rhng '
    ratenames(irdgn)   = 'rdgn '
    ratenames(irdpg)   = 'rdpg '
    ratenames(irhegp)  = 'rhegp'

    ! for weak reactions
    ratenames(irpen)   = 'rpen '
    ratenames(irnep)   = 'rnep '
    ratenames(irn56ec) = 'r56ec'

    ! ppchain
    ratenames(irpp)    = 'rpp  '
    ratenames(ir33)    = 'r33  '
    ratenames(irhe3ag) = 'rhe3ag'

    ! cno cycles
    ratenames(ircpg)   = 'rcpg '
    ratenames(irnpg)   = 'rnpg '
    ratenames(iropg)   = 'ropg '
    ratenames(ifa)     = 'rfa  '
    ratenames(ifg)     = 'rfg  '
    ratenames(irnag)   = 'rnag '


    ! the dummy links
    ratenames(irr1)   = 'r1   '
    ratenames(irs1)   = 's1   '
    ratenames(irt1)   = 't1   '
    ratenames(iru1)   = 'u1   '
    ratenames(irv1)   = 'v1   '
    ratenames(irw1)   = 'w1   '
    ratenames(irx1)   = 'x1   '

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

  end subroutine actual_network_init


  subroutine actual_network_finalize

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    if (allocated(aion)) then
       deallocate(aion)
    endif
    if (allocated(zion)) then
       deallocate(zion)
    endif
    if (allocated(nion)) then
       deallocate(nion)
    endif
    if (allocated(bion)) then
       deallocate(bion)
    endif
    if (allocated(mion)) then
       deallocate(mion)
    endif
    if (allocated(wion)) then
       deallocate(wion)
    endif

  end subroutine actual_network_finalize

end module actual_network
