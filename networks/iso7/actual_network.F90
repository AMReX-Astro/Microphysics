module actual_network

  use network_properties
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: ihe4  = 1
  integer, parameter :: ic12  = 2
  integer, parameter :: io16  = 3
  integer, parameter :: ine20 = 4
  integer, parameter :: img24 = 5
  integer, parameter :: isi28 = 6
  integer, parameter :: ini56 = 7

  real(rt)        , allocatable :: bion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: bion, mion, wion
#endif

  character (len=32), parameter :: network_name = "iso7"

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

  integer, parameter :: nrates  = 17
  integer, parameter :: num_rate_groups = 4

  integer, parameter :: ircag   = 1
  integer, parameter :: iroga   = 2
  integer, parameter :: ir3a    = 3
  integer, parameter :: irg3a   = 4
  integer, parameter :: ir1212  = 5
  integer, parameter :: ir1216  = 6
  integer, parameter :: ir1616  = 7
  integer, parameter :: iroag   = 8
  integer, parameter :: irnega  = 9
  integer, parameter :: irneag  = 10
  integer, parameter :: irmgga  = 11
  integer, parameter :: irmgag  = 12
  integer, parameter :: irsiga  = 13
  integer, parameter :: ircaag  = 14
  integer, parameter :: irtiga  = 15
  integer, parameter :: irsi2ni = 16
  integer, parameter :: irni2si = 17

  character (len=20), save :: ratenames(nrates)

contains

  subroutine actual_network_init

    implicit none

    call network_properties_init()

    allocate(bion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))


    ! Set the binding energy of the element
    bion(ihe4)  = 28.29603e0_rt
    bion(ic12)  = 92.16294e0_rt
    bion(io16)  = 127.62093e0_rt
    bion(ine20) = 160.64788e0_rt
    bion(img24) = 198.25790e0_rt
    bion(isi28) = 236.53790e0_rt
    bion(ini56) = 484.00300e0_rt

    ! Set the mass
    mion(:) = nion(:) * mn + zion(:) * (mp + me) - bion(:) * mev2gr

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    wion(:) = aion(:)

    ! set the names of the reaction rates
    ratenames(ircag)   = 'rcag '
    ratenames(iroga)   = 'roga '
    ratenames(ir3a)    = 'r3a  '
    ratenames(irg3a)   = 'rg3a '    ! inverse rate
    ratenames(ir1212)  = 'r1212'
    ratenames(ir1216)  = 'r1216'
    ratenames(ir1616)  = 'r1616'
    ratenames(iroag)   = 'roag '
    ratenames(irnega)  = 'rnega'
    ratenames(irneag)  = 'rneag'
    ratenames(irmgga)  = 'rmgga'
    ratenames(irmgag)  = 'rmgag'
    ratenames(irsiga)  = 'rsiga'
    ratenames(ircaag)  = 'rcaag'
    ratenames(irtiga)  = 'rtiga'

    ratenames(irsi2ni) = 'rsi2ni'
    ratenames(irni2si) = 'rni2si'

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    call network_properties_init()

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
