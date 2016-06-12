! Setup a 3D grid of smoothly varying rho, T, and user-defined X.  Then
! call react_state() on the grid and output the results.

program test_react

  use BoxLib
  use bl_constants_module
  use bl_types
  use bl_space
  use f2kcli
  use box_util_module
  use ml_layout_module
  use multifab_module
  use variables
  use probin_module, only: dens_min, dens_max, &
                           temp_min, temp_max, test_set, dt
  use runtime_init_module
  use eos_module
  use eos_type_module
  use burn_type_module
  use actual_burner_module
  use microphysics_module
  use network
  use util_module
  use variables

  !Local variables
  implicit none

  ! Conventional fluid state multifabs
  type(multifab) , allocatable :: s(:)
  
  real(kind=dp_t) :: dx(1, MAX_SPACEDIM)

  logical :: pmask(MAX_SPACEDIM)

  type(ml_layout) :: mla
  type(ml_boxarray) :: mba

  integer :: i, n
  integer :: ii, jj, kk
  integer :: nrho, nT, nX

  integer :: dm, nlevs

  type(plot_t) :: pf

  real(kind=dp_t), pointer :: sp(:,:,:,:)                                  
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)                     
  integer :: domlo(MAX_SPACEDIM), domhi(MAX_SPACEDIM)
  
  type (burn_t) :: burn_state_in, burn_state_out
  type (eos_t) :: eos_state

  real (kind=dp_t) :: dens_zone, temp_zone
  real (kind=dp_t) :: dlogrho, dlogT
  real (kind=dp_t), allocatable :: xn_zone(:, :)

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)                                          


  call runtime_init()
  call init_variables(pf)

  ! initialize a grid -- since we are not doing anything hydroy, set the
  ! periodic mask to true
  pmask(:) = .true.
  call read_a_hgproj_grid(mba, test_set)
  call ml_layout_build(mla, mba, pmask)

  nlevs = mla % nlevel
  if (nlevs /= 1) then
      call bl_error("ERROR: only 1 level of refinement currently supported")
   endif

  dm = mla % dim
  if (dm /= 3) then
     call bl_error("ERROR: we require dm = 3")
  endif

  ! we don't care about dx -- we have no physical size
  dx(1,:) = ONE

  ! microphysics
  call microphysics_init()

  ! we'll store everything in a multifab -- inputs and outputs
  call init_variables(pf)

  allocate(s(nlevs))

  do n = 1,nlevs
    call multifab_build(s(n), mla%la(n), pf % n_plot_comps, 0)
  end do

  nrho = extent(mla%mba%pd(1),1)
  nT = extent(mla%mba%pd(1),2)
  nX = extent(mla%mba%pd(1),3)

  dlogrho = (log10(dens_max) - log10(dens_min))/(nrho - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(nT - 1)

  ! read from the input file to get all the species data
  domlo = lwb(get_pd(get_layout(s(1))))
  domhi = upb(get_pd(get_layout(s(1))))

  allocate(xn_zone(nspec, 0:nX-1))   ! this assumes that lo(3) = 0

  call get_xn(xn_zone, domlo(3), domhi(3))

  n = 1  ! single level assumption
  do i = 1, nfabs(s(n))
     sp => dataptr(s(n), i)

     lo = lwb(get_box(s(n), i))
     hi = upb(get_box(s(n), i))

     do kk = lo(3), hi(3)
        do jj = lo(2), hi(2)
           temp_zone = 10.0_dp_t**(log10(temp_min) + dble(jj)*dlogT)

           do ii = lo(1), hi(1)
              dens_zone = 10.0_dp_t**(log10(dens_min) + dble(ii)*dlogrho)

              burn_state_in % rho = dens_zone
              burn_state_in % T = temp_zone

              burn_state_in % xn(:) = xn_zone(:, kk)

              call burn_to_eos(burn_state_in, eos_state)
              call eos(eos_input_rt, eos_state)
              call eos_to_burn(eos_state, burn_state_in)

              call actual_burner(burn_state_in, burn_state_out, dt, ZERO)

              ! store
              sp(ii, jj, kk, pf % irho) = dens_zone
              sp(ii, jj, kk, pf % itemp) = temp_zone
              sp(ii, jj, kk, pf % ispec: pf % ispec-1+nspec) = burn_state_out % xn(:)
              sp(ii, jj, kk, pf % irodot: pf % irodot-1+nspec) = &
                   (burn_state_out % xn(:) - burn_state_in % xn(:)) / dt
              sp(ii, jj, kk, pf % irho_hnuc) = &
                   dens_zone * (burn_state_out % e - burn_state_in % e) / dt

           enddo
        enddo
     enddo

  enddo


  ! output


  ! if you (or a subroutine) built it, destroy it!
  do n = 1,nlevs
    call destroy(s(n))
  end do

  call destroy(mla)

  deallocate(s)

  call runtime_close()

  ! end boxlib
  call bl_prof_glean("bl_prof_res")                                             
  call bl_prof_finalize()                                                       
  call boxlib_finalize()                                                        

end program test_react
