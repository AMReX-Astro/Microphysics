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
                           temp_min, temp_max, test_set, &
                           metalicity_max, &
                           small_temp, small_dens
  use runtime_init_module
  use microphysics_module
  use network
  use eos_type_module
  use eos_module
  use variables
  use fabio_module
  use build_info_module

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
  
  type (eos_t) :: eos_state, eos_state_reference

  real(kind=dp_t) :: temp_zone, dens_zone, metalicity
  real(kind=dp_t) :: dlogrho, dlogT, dmetal
  real(kind=dp_t) :: xn_zone(nspec)

  character (len=128) :: out_name

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
  dmetal    = (metalicity_max  - ZERO)/(nX - 1)

  ! read from the input file to get all the species data
  domlo = lwb(get_pd(get_layout(s(1))))
  domhi = upb(get_pd(get_layout(s(1))))

  n = 1  ! single level assumption
  do i = 1, nfabs(s(n))
     sp => dataptr(s(n), i)

     lo = lwb(get_box(s(n), i))
     hi = upb(get_box(s(n), i))

     !$OMP PARALLEL DO PRIVATE(ii,jj,kk,metalicity,temp_zone,dens_zone,eos_state_reference,xn_zone) &
     !$OMP FIRSTPRIVATE (eos_state)
     do kk = lo(3), hi(3)
        ! set the composition -- approximately solar
        metalicity = ZERO + dble(kk)*dmetal
        xn_zone(:) = metalicity/(nspec - 2)   ! all but H, He
        xn_zone(ih1)  = 0.75_dp_t - HALF*metalicity
        xn_zone(ihe4) = 0.25_dp_t - HALF*metalicity

        do jj = lo(2), hi(2)
           temp_zone = 10.0_dp_t**(log10(temp_min) + dble(jj)*dlogT)

           do ii = lo(1), hi(1)
              dens_zone = 10.0_dp_t**(log10(dens_min) + dble(ii)*dlogrho)

              eos_state % rho = dens_zone
              eos_state % T = temp_zone
              eos_state % xn(:) = xn_zone(:)

              ! store default state
              sp(ii, jj, kk, pf % irho) = dens_zone
              sp(ii, jj, kk, pf % itemp) = temp_zone
              sp(ii, jj, kk, pf % ispec: pf % ispec-1+nspec) = xn_zone(:)


              ! call EOS using rho, T
              call eos(eos_input_rt, eos_state)

              eos_state_reference = eos_state

              sp(ii, jj, kk, pf % ih) = eos_state % h
              sp(ii, jj, kk, pf % ie) = eos_state % e
              sp(ii, jj, kk, pf % ip) = eos_state % p
              sp(ii, jj, kk, pf % is) = eos_state % s


              ! call EOS using rho, h

              ! reset T to give it some work to do
              eos_state % T = 100.d0

              call eos(eos_input_rh, eos_state)

              sp(ii, jj, kk, pf % ierr_T_eos_rh) = &
                   abs(eos_state % T - temp_zone)/temp_zone

              eos_state = eos_state_reference


              ! call EOS using T, p

              ! reset rho to give it some work to do
              eos_state % rho = 1.d0

              call eos(eos_input_tp, eos_state)

              sp(ii, jj, kk, pf % ierr_rho_eos_tp) = &
                   abs(eos_state % rho - dens_zone)/dens_zone

              eos_state = eos_state_reference


              ! call EOS using r, p

              ! reset T to give it some work to do
              eos_state % T = 100.d0

              call eos(eos_input_rp, eos_state)

              sp(ii, jj, kk, pf % ierr_T_eos_rp) = &
                   abs(eos_state % T - temp_zone)/temp_zone

              eos_state = eos_state_reference



              ! call EOS using r, e

              ! reset T to give it some work to do
              eos_state % T = 100.d0

              call eos(eos_input_re, eos_state)

              sp(ii, jj, kk, pf % ierr_T_eos_re) = &
                   abs(eos_state % T - temp_zone)/temp_zone

              eos_state = eos_state_reference



              ! call EOS using p, s

              ! reset T and rho to give it some work to do
              eos_state % T = 100.d0
              eos_state % rho = 1.d0


              ! some EOSes don't have physically valid treatments
              ! of entropy throughout the entire rho-T plane
              if (eos_state%s > ZERO) then

                 call eos(eos_input_ps, eos_state)
                 
                 ! store the thermodynamic state
                 sp(ii, jj, kk, pf % ierr_T_eos_ps) = &
                      abs(eos_state % T - temp_zone)/temp_zone
                 sp(ii, jj, kk, pf % ierr_rho_eos_ps) = &
                      abs(eos_state % rho - dens_zone)/dens_zone
                 
              else
                 sp(ii, jj, kk, pf % ierr_T_eos_ps) = ZERO
                 sp(ii, jj, kk, pf % ierr_rho_eos_ps) = ZERO

              endif

              eos_state = eos_state_reference

              
              ! call EOS using p, h

              ! reset T and rho to give it some work to do
              eos_state % T = 100.d0
              eos_state % rho = 1.d0

              call eos(eos_input_ph, eos_state)
                 
              sp(ii, jj, kk, pf % ierr_T_eos_ph) = &
                   abs(eos_state % T - temp_zone)/temp_zone
              sp(ii, jj, kk, pf % ierr_rho_eos_ph) = &
                   abs(eos_state % rho - dens_zone)/dens_zone
                 
              eos_state = eos_state_reference


              ! call EOS using T, h
              ! this doesn't work for all EOSes (where h doesn't depend on T)

              if (.not. trim(eos_dir) == "gamma_law_general") then
                 ! reset rho to give it some work to do -- for helmeos, h is not
                 ! monotonic, so we only perturb rho slightly here
                 eos_state % rho = 0.9 * eos_state % rho

                 call eos(eos_input_th, eos_state)
                 
                 sp(ii, jj, kk, pf % ierr_rho_eos_th) = &
                      abs(eos_state % rho - dens_zone)/dens_zone
                 
                 eos_state = eos_state_reference

              else
                 sp(ii, jj, kk, pf % ierr_rho_eos_th) = ZERO
              endif

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

  enddo


  ! output
  out_name = "test_eos." // trim(eos_dir)

  call fabio_ml_multifab_write_d(s, mla%mba%rr(:,1), out_name, names=pf%names)
  call write_job_info(out_name, mla%mba)


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
