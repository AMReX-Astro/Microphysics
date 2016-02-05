module actual_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network
  use actual_burner_data
  use burn_type_module

  implicit none

  ! Conversion factor for the nuclear energy generation rate.

  double precision, parameter, private :: avo = 6.0221417930d23
  double precision, parameter, private :: c_light = 2.99792458d10
  double precision, parameter, private :: enuc_conv2 = -avo*c_light*c_light

contains

  subroutine actual_burner(state_in, state_out, dt, time)

    use integration_module, only: do_burn

    implicit none

    type (burn_t),       intent(in   ) :: state_in
    type (burn_t),       intent(inout) :: state_out
    double precision,    intent(in   ) :: dt, time

    call do_burn(state_in, state_out, dt, time)

  end subroutine actual_burner



  subroutine actual_burner_init()

    use integration_module, only: integration_init
    use rates_module, only: rates_init
    use screening_module, only: screening_init

    implicit none

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

    call integration_init()

    call rates_init()

    call set_up_screening_factors()

    call screening_init()

  end subroutine actual_burner_init



  ! Compute and store the more expensive screening factors  

  subroutine set_up_screening_factors()

    use screening_module, only: add_screening_factor
    use network, only: aion, zion

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
