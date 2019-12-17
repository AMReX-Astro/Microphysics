module actual_network

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt)        , parameter :: MeV2erg = 1.60217646e-6_rt
  real(rt)        , parameter :: N_A = 6.0221415e23_rt

  integer, parameter :: nspec = 10
  integer, parameter :: nspec_evolve = 10
  integer, parameter :: naux  = 0

  integer, parameter :: ic12  = 1
  integer, parameter :: io14  = 2
  integer, parameter :: io15  = 3
  integer, parameter :: io16  = 4
  integer, parameter :: if17  = 5
  integer, parameter :: img22 = 6
  integer, parameter :: is30  = 7
  integer, parameter :: ini56 = 8
  integer, parameter :: ihe4  = 9
  integer, parameter :: ih1   = 10

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(rt)        , allocatable :: aion(:), zion(:), ebin(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, ebin
#endif

  character (len=32), parameter :: network_name = "rprox"

  ! Rates data

  integer, parameter :: nrates    = 18
  integer, parameter :: num_rate_groups = 3

  integer, parameter :: irlambCNO = 1
  integer, parameter :: irag15o   = 2
  integer, parameter :: irr1      = 3
  integer, parameter :: irag16o   = 4
  integer, parameter :: irpg16o   = 5
  integer, parameter :: irpg17f   = 6
  integer, parameter :: irgp17f   = 7
  integer, parameter :: irlambda2 = 8
  integer, parameter :: irap14o   = 9
  integer, parameter :: irs1      = 10
  integer, parameter :: irlambda1 = 11
  integer, parameter :: ir3a      = 12
  integer, parameter :: irpg12c   = 13
  integer, parameter :: irwk14o   = 14
  integer, parameter :: irwk17f   = 15
  integer, parameter :: irwk15o   = 16
  integer, parameter :: irLweak   = 17
  integer, parameter :: irla2     = 18

  character (len=10), save :: reac_names(nrates)

  integer, parameter :: dlambCNOdh1   = 1
  integer, parameter :: drs1dhe4      = 2
  integer, parameter :: drr1dh1       = 3
  integer, parameter :: dlambda1dhe4  = 4
  integer, parameter :: dlambda2dhe4  = 5
  integer, parameter :: delta1        = 6
  integer, parameter :: delta2        = 7
  integer, parameter :: r56eff        = 8
  integer, parameter :: dr56effdt     = 9

contains

  subroutine actual_network_init

    use amrex_constants_module

    ! set the names
    spec_names(ic12)  = "carbon-12"
    spec_names(io14)  = "oxygen-14"
    spec_names(io15)  = "oxygen-15"
    spec_names(io16)  = "oxygen-16"
    spec_names(if17)  = "flourine-17"
    spec_names(img22) = "magnesium-22"
    spec_names(is30)  = "sulfur-30"
    spec_names(ini56) = "nickel-56"
    spec_names(ihe4)  = "helium-4"
    spec_names(ih1)   = "hydrogen-1"

    short_spec_names(ic12)  = "C12"
    short_spec_names(io14)  = "O14"
    short_spec_names(io15)  = "O15"
    short_spec_names(io16)  = "O16"
    short_spec_names(if17)  = "F17"
    short_spec_names(img22) = "Mg22"
    short_spec_names(is30)  = "S30"
    short_spec_names(ini56) = "Ni56"
    short_spec_names(ihe4)  = "He4"
    short_spec_names(ih1)   = "H1"    

    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(ebin(nspec))

    ! set the species properties
    aion(ic12)  = TWELVE
    aion(io14)  = 14.0_rt
    aion(io15)  = 15.0_rt
    aion(io16)  = 16.0_rt
    aion(if17)  = 17.0_rt
    aion(img22) = 22.0_rt
    aion(is30)  = 30.0_rt
    aion(ini56) = 56.0_rt
    aion(ihe4)  = FOUR
    aion(ih1)   = ONE

    zion(ic12)  = SIX
    zion(io14)  = EIGHT
    zion(io15)  = EIGHT
    zion(io16)  = EIGHT
    zion(if17)  = NINE
    zion(img22) = TWELVE
    zion(is30)  = 16.0_rt
    zion(ini56) = 28.0_rt
    zion(ihe4)  = TWO
    zion(ih1)   = ONE

    ! Our convention is that binding energy is negative.  The
    ! following are the binding energies in MeV.
    ebin(ic12)  = 92.16279_rt
    ebin(io14)  = 98.7325_rt
    ebin(io15)  = 111.9569_rt
    ebin(io16)  = 127.6207_rt
    ebin(if17)  = 128.2211_rt
    ebin(img22) = 168.5768_rt
    ebin(is30)  = 243.6866_rt
    ebin(ini56) = 483.995_rt
    ebin(ihe4)  = 28.29599_rt
    ebin(ih1)   = ZERO

    ! convert to erg / g by multiplying by N_A / aion and converting to erg
    ebin = -ebin * N_A * MeV2erg / aion

    reac_names(irlambCNO) = "rlambdaCNO"
    reac_names(irag15o)   = "rag15o"
    reac_names(irr1)      = "rr1"
    reac_names(irag16o)   = "rag16o"
    reac_names(irpg16o)   = "rpg16o"
    reac_names(irpg17f)   = "rpg17f"
    reac_names(irgp17f)   = "rgp17f"
    reac_names(irlambda2) = "rlambda2"
    reac_names(irap14o)   = "rap14o"
    reac_names(irs1)      = "rs1"
    reac_names(irlambda1) = "rlambda1"
    reac_names(ir3a)      = "r3a"
    reac_names(irpg12c)   = "rpg12c"
    reac_names(irwk14o)   = "wk14o"
    reac_names(irwk17f)   = "wk17f"
    reac_names(irwk15o)   = "wk15o"
    reac_names(irLweak)   = "Lweak"
    reac_names(irla2)     = "la2"

  end subroutine actual_network_init



  subroutine actual_network_finalize

    implicit none

    if (allocated(aion)) then
       deallocate(aion)
    endif
    if (allocated(zion)) then
       deallocate(zion)
    endif
    if (allocated(ebin)) then
       deallocate(ebin)
    endif

  end subroutine actual_network_finalize

end module actual_network
