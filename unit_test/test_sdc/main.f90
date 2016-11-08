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
  use variables, only: init_variables, finalize_variables, plot_t
  use probin_module, only: dens_min, dens_max, &
                           temp_min, temp_max, test_set, tmax, run_prefix, &
                           small_temp, small_dens, do_acc
  use runtime_init_module
  use sdc_type_module
  use microphysics_module
  use integrator_module, only: integrator
  use eos_type_module, only : eos_t, eos_get_small_temp, eos_get_small_dens
  use eos_module, only: eos, eos_input_rt
  use network, only: nspec
  use util_module
  use fabio_module
  use build_info_module
  use parallel, only : parallel_wtime

  implicit none

  ! Conventional fluid state multifabs
  type(multifab) , allocatable :: s(:)

  real(kind=dp_t) :: dx(1, MAX_SPACEDIM)

  logical :: pmask(MAX_SPACEDIM)

  type(ml_layout) :: mla
  type(ml_boxarray) :: mba

  integer :: i, j, n
  integer :: ii, jj, kk
  integer :: nrho, nT, nX

  integer :: dm, nlevs

  integer :: n_rhs_min, n_rhs_max, n_rhs_avg

  type(plot_t) :: pf

  integer :: itemp, irho, ispec, ispec_old, irodot, irho_hnuc

  real(kind=dp_t), pointer :: sp(:,:,:,:)

  real(kind=dp_t), allocatable :: state(:,:,:,:)

  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: domlo(MAX_SPACEDIM), domhi(MAX_SPACEDIM)

  type (sdc_t) :: sdc_state_in, sdc_state_out
  type (eos_t) :: eos_state

  real (kind=dp_t) :: dlogrho, dlogT
  real (kind=dp_t), allocatable :: xn_zone(:, :)

  real (kind=dp_t) :: sum_X

  real (kind=dp_t) :: start_time, end_time

  character (len=256) :: out_name

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)


  call runtime_init(.true.)

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
  call microphysics_init(small_temp=small_temp, small_dens=small_dens)

  call eos_get_small_temp(small_temp)
  print *, "small_temp = ", small_temp

  call eos_get_small_dens(small_dens)
  print *, "small_dens = ", small_dens

  ! we'll store everything in a multifab -- inputs and outputs
  call init_variables(pf)

  allocate(s(nlevs))

  do n = 1,nlevs
    call multifab_build(s(n), mla%la(n), pf % n_plot_comps, 0)
  end do

  nrho = extent(mla%mba%pd(1),1)
  nT = extent(mla%mba%pd(1),2)
  nX = extent(mla%mba%pd(1),3)

  allocate(state(0:nrho-1, 0:nT-1, 0:nX-1, pf % n_plot_comps))

  dlogrho = (log10(dens_max) - log10(dens_min))/(nrho - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(nT - 1)

  ! read from the input file to get all the species data
  domlo = lwb(get_pd(get_layout(s(1))))
  domhi = upb(get_pd(get_layout(s(1))))

  allocate(xn_zone(nspec, 0:nX-1))   ! this assumes that lo(3) = 0

  call get_xn(xn_zone, domlo(3), domhi(3))

  ! normalize -- just in case
  do kk = domlo(3), domhi(3)
     sum_X = sum(xn_zone(:, kk))
     xn_zone(:, kk) = xn_zone(:, kk)/sum_X
  enddo

  ! GPU doesn't like derived-types with bound procedures
  itemp = pf % itemp
  irho = pf % irho
  ispec = pf % ispec
  ispec_old = pf % ispec_old
  irodot = pf % irodot
  irho_hnuc = pf % irho_hnuc

  n = 1  ! single level assumption

  n_rhs_avg = 0
  n_rhs_max = -100000000
  n_rhs_min = 100000000

  do i = 1, nfabs(s(n))
     sp => dataptr(s(n), i)

     lo = lwb(get_box(s(n), i))
     hi = upb(get_box(s(n), i))

     ! First, construct the input state in a separate loop.

     do kk = lo(3), hi(3)
        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)

              state(ii, jj, kk, pf % itemp) = 10.0_dp_t**(log10(temp_min) + dble(jj)*dlogT)
              state(ii, jj, kk, pf % irho) = 10.0_dp_t**(log10(dens_min) + dble(ii)*dlogrho)
              state(ii, jj, kk, pf%ispec_old:pf%ispec_old+nspec-1) = max(xn_zone(:, kk), 1.e-10_dp_t)

           enddo
        enddo
     enddo

     ! Set up a timer for the burn.

     start_time = parallel_wtime()

     !$OMP PARALLEL DO PRIVATE(ii,jj,kk,j) &
     !$OMP PRIVATE(burn_state_in, burn_state_out) &
     !$OMP REDUCTION(+:n_rhs_avg) REDUCTION(MAX:n_rhs_max) REDUCTION(MIN:n_rhs_min) &
     !$OMP SCHEDULE(DYNAMIC,1)

     !$acc data copyin(temp_min, dlogT, dens_min, dlogrho, lo, hi, tmax) &
     !$acc      copyin(itemp, irho, ispec, ispec_old, irodot, irho_hnuc) &
     !$acc      copy(state(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)) if (do_acc == 1)

     !$acc parallel reduction(+:n_rhs_avg) reduction(max:n_rhs_max) reduction(min:n_rhs_min) if (do_acc == 1)

     !$acc loop gang vector collapse(3) &
     !$acc private(burn_state_in, burn_state_out, ii, jj, kk, j)

     do kk = lo(3), hi(3)
        do jj = lo(2), hi(2)
           do ii = lo(1), hi(1)

              ! populate the SDC state.  Our strategy for the unit test is to choose
              ! the advective terms to all be zero and to choose the velocity to be
              ! Mach = 0.1

              ! call the EOS first to get the sound speed
              eos_state % rho = state(ii, jj, kk, irho)
              eos_state % T = state(ii, jj, kk, itemp)
              eos_state % xn(:) = state(ii, jj, kk, ispec_old:ispec_old-1+nspec)

              call eos(eos_input_rt, eos_state)

              sdc_state_in % y(SRHO) = state(ii, jj, kk, irho)

              ! we will pick velocities to be 10% of the sound speed
              sdc_state_in % y(SMX:SMZ) = sdc_state_in % y(SRHO) * 0.1 * eos_state % cs

              sdc_state_in % y(SEINT) = sdc_state_in % y(SRHO) * eos_state % e
              sdc_state_in % y(SEDEN) = sdc_state_in % y(SEINT) + &
                   HALF*sum(sdc_state_in % y(SMX:SMZ)**2)/sdc_state_in % y(SRHO)
              sdc_state_in % y(SFS:SFS-1+nspec) = sdc_state_in % y(SRHO) * eos_state % xn(:)

              ! zero out the advective terms
              sdc_state_in % ydot_a(:) = ZERO

              ! need to set T_from_eden

              call integrator(sdc_state_in, sdc_state_out, tmax, ZERO)


              n_rhs_avg = n_rhs_avg + sdc_state_out % n_rhs
              n_rhs_min = min(n_rhs_min, sdc_state_out % n_rhs)
              n_rhs_max = max(n_rhs_max, sdc_state_out % n_rhs)

           enddo
        enddo
     enddo
     !$acc end parallel
     !$acc end data

     !$OMP END PARALLEL DO

     ! End the timer and print the results.

     end_time = parallel_wtime()

     print *, "Execution time: ", end_time - start_time

     sp(:,:,:,:) = state(:,:,:,:)
  enddo

  ! note: integer division
  n_rhs_avg = n_rhs_avg/(nT*nrho*nX)

  print *, "RHS stats:"
  print *, "  min: ", n_rhs_min
  print *, "  avg: ", n_rhs_avg
  print *, "  max: ", n_rhs_max

  ! output
  out_name = trim(run_prefix) // "test_react." // trim(integrator_dir)

  call fabio_ml_multifab_write_d(s, mla%mba%rr(:,1), trim(out_name), names=pf%names)

  call write_job_info(out_name, mla%mba)


  ! if you (or a subroutine) built it, destroy it!
  do n = 1,nlevs
    call destroy(s(n))
  end do

  call destroy(mla)

  call finalize_variables(pf)

  deallocate(state)
  deallocate(s)
  deallocate(xn_zone)

  call runtime_close()


  call microphysics_finalize()

  ! end boxlib
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

end program test_react
