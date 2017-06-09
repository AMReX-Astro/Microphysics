! Setup a 3D grid of smoothly varying rho, T, and user-defined X.  Then
! call react_state() on the grid and output the results.

program test_react
#ifdef CUDA
  use cudafor
#endif
  use react_zones_module, only: pfidx_t, react_zones
  
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
                           temp_min, temp_max, test_set, run_prefix, &
                           small_temp, small_dens, do_acc
  use runtime_init_module
  use burn_type_module
  use actual_burner_module, only : actual_burner
  use microphysics_module
  use eos_type_module, only : mintemp, mindens
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

  real(kind=dp_t), &
#ifdef CUDA
       managed, &
#endif
       allocatable :: state(:,:,:,:)

  integer, &
#ifdef CUDA
       managed, &
#endif
       allocatable :: lo(:), hi(:)
  
  integer :: domlo(MAX_SPACEDIM), domhi(MAX_SPACEDIM)

  real (kind=dp_t) :: dlogrho, dlogT
  real (kind=dp_t), allocatable :: xn_zone(:, :)

  real (kind=dp_t) :: sum_X

  real (kind=dp_t) :: start_time, end_time

  character (len=256) :: out_name

  type(pfidx_t), &
#ifdef CUDA
       managed, &
#endif
       allocatable :: pfidx
  
#ifdef CUDA  
  type(dim3)    :: cuGrid, cuThreadBlock
#endif
  
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
  call microphysics_init(small_temp=mintemp, small_dens=mindens)

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
  allocate(pfidx)
  pfidx % itemp = pf % itemp
  pfidx % irho = pf % irho
  pfidx % ispec = pf % ispec
  pfidx % ispec_old = pf % ispec_old
  pfidx % irodot = pf % irodot
  pfidx % irho_hnuc = pf % irho_hnuc
  pfidx % endrho = nrho - 1
  pfidx % endT = nT - 1
  pfidx % endX = nX - 1
  pfidx % ncomps = pf % n_plot_comps

  n = 1  ! single level assumption

  n_rhs_avg = 0
  n_rhs_max = -100000000
  n_rhs_min = 100000000

  ! Allocate lo, hi
  allocate(lo(MAX_SPACEDIM))
  allocate(hi(MAX_SPACEDIM))
  
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

     write(*,*) 'lo = ', lo
     write(*,*) 'hi = ', hi
     
#ifdef CUDA
     ! Set up CUDA parameters
     cuThreadBlock = dim3(32, 8, 1)
     cuGrid = dim3(&
          ceiling(real(hi(1)-lo(1)+1)/cuThreadBlock%x), &
          ceiling(real(hi(2)-lo(2)+1)/cuThreadBlock%y), &
          ceiling(real(hi(3)-lo(3)+1)/cuThreadBlock%z))
     
     write(*,*) 'cuGrid % x = ', cuGrid % x
     write(*,*) 'cuGrid % y = ', cuGrid % y
     write(*,*) 'cuGrid % z = ', cuGrid % z
     
     ! React the zones using CUDA     
     call react_zones<<<cuGrid,cuThreadBlock>>>(state, pfidx, lo, hi)
#else
     call react_zones(state, pfidx, lo, hi)
#endif

     !! Do reduction on statistics
     ! n_rhs_avg = n_rhs_avg + burn_state_out % n_rhs
     ! n_rhs_min = min(n_rhs_min, burn_state_out % n_rhs)
     ! n_rhs_max = max(n_rhs_max, burn_state_out % n_rhs)
     
     ! End the timer and print the results.

     end_time = parallel_wtime()

     print *, "Execution time: ", end_time - start_time

     sp(:,:,:,:) = state(:,:,:,:)
  enddo

  ! note: integer division
  ! n_rhs_avg = n_rhs_avg/(nT*nrho*nX)

  ! print *, "RHS stats:"
  ! print *, "  min: ", n_rhs_min
  ! print *, "  avg: ", n_rhs_avg
  ! print *, "  max: ", n_rhs_max

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
  
  deallocate(pfidx)
  deallocate(state)
  deallocate(s)
  deallocate(xn_zone)

  call runtime_close()


  call microphysics_finalize()

  ! end boxlib
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

  call send_success_return_code()
  
end program test_react
