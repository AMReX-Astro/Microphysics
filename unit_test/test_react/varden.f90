! Setup a 3D grid of smoothly varying rho, T, and user-defined X.  Then
! call react_state() on the grid and output the results.

subroutine varden()

  use BoxLib
  use f2kcli
  use ml_layout_module
  use multifab_module
  use variables
  use probin_module, only: small_temp, &
                           dens_min, base_cutoff_density
  use runtime_init_module
  use bl_constants_module
  use bl_types
  use network
  use eos_module
  use varden_aux
  use variables

  !Local variables
  implicit none

  ! Conventional fluid state multifabs
  type(multifab) , allocatable :: s(:), snew(:)

  ! react_state output
  type(multifab) , allocatable :: rho_omegadot(:), rho_Hnuc(:), rho_Hext(:)

  real(kind=dp_t), allocatable :: tempbar(:,:), pbar(:,:)

  real(kind=dp_t), pointer :: dx(:,:)

  type(ml_layout) :: mla
  type(bc_tower)  :: bct

  character(len=100) :: temp_buf

  real(kind=dp_t) :: m_in, m_out

  integer :: i, n, res
  integer :: ii, jj, kk
  integer :: nlevs

  logical :: dbo, dho

  type(plot_t) :: pf

  call runtime_init()
  call init_variables()

  ! initialize a grid -- since we are not doing anything hydroy, set the
  ! periodic mask to true
  pmask(:) = .true.
  call read_a_hgproj_grid(mba, test_set)
  call ml_layout_build(mla, mba, pmask)

  nlevs = mla % nlevel
  if (nlevs /= 1) then
      call bl_error("ERROR: only 1 level of refinement currently supported")
   endif

  dm = mla % dm
  if (dm /= 3) then
     call bl_error("ERROR: we require dm = 3")
  endif

  ! we don't care about dx -- we have no physical size
  dx(1,:) = ONE

  ! microphysics
  call network_init()
  call eos_init()

  ! we'll store everything in a multifab -- inputs and outputs
  call init_variables(pf)

  allocate(s(nlevs))

  do n = 1,nlevs
    call multifab_build(s(n), mla%la(n), pf % n_plot_comps, 0)
  end do

  nrho = (extent(mla%mba%pd(1),1)
  nT = (extent(mla%mba%pd(1),2)
  nX = (extent(mla%mba%pd(1),3)

  dlogrho = (log10(dens_max) - log10(dens_min))/(nrho - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(nT - 1)

  ! read from the input file to get all the species data
  domlo = lwb(get_pd(get_layout(s(1))))
  domhi = upb(get_pd(get_layout(s(1))))
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

              burn_state % rho = dens_zone
              burn_state % T = temp_zone

              burn_state % xn(:) = xn_zone(:, kk)

              call burn_to_eos(burn_state_in, eos_state)
              call eos(eos_input_rt, eos_state)
              call eos_to_burn(eos_state, burn_state_in)

              call actual_burner(burn_state_in, burn_state_out, dt, ZERO)

              ! store
              sp(ii, jj, kk, irho) = dens_zone
              sp(ii, jj, kk, itemp) = temp_zone
              sp(ii, jj, kk, ispec:ispec-1+nspec) = burn_state_out % xn(:)
              sp(ii, jj, kk, irodot:irodot-1+nspec) = &
                   (burn_state_out % xn(:) - burn_state_in % xn(:)) / dt
              sp(ii, jj, kk, irho_hnuc) = &
                   dens_zone * (burn_state_out % e - burn_state_in % e) / dt

           enddo
        enddo
     enddo

  enddo


  ! output

  call varden_close()

  ! if you (or a subroutine) built it, destroy it!
  do n = 1,nlevs
    call destroy(s(n))
  end do

  call destroy(mla)

  deallocate(s)

  call runtime_close()

end subroutine varden
