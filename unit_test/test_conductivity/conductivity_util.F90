subroutine do_conductivity(lo, hi, &
                           sp, s_lo, s_hi, npts) bind(C, name="do_conductivity")

  use variables
  use network
  use eos_type_module
  use eos_module
  use microphysics_type_module
  use extern_probin_module
  use actual_conductivity_module, only: actual_conductivity

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: s_lo(3), s_hi(3)
  real(rt), intent(inout) :: sp(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), p % n_plot_comps)
  integer, intent(in), value :: npts

  real(rt) :: dlogrho, dlogT, dmetal
  real(rt) :: metalicity, temp_zone, dens_zone
  real(rt) :: xn_zone(nspec)

  type(eos_t) :: eos_state
  type(eos_t) :: eos_state_reference

  integer :: ii, jj, kk

  !$gpu
  
  dlogrho = (log10(dens_max) - log10(dens_min))/(npts - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(npts - 1)
  dmetal    = (metalicity_max  - ZERO)/(npts - 1)

  do kk = lo(3), hi(3)

     ! set the composition -- approximately solar
     metalicity = ZERO + dble(kk)*dmetal
     xn_zone(:) = metalicity/(nspec - 2)   ! all but H, He
     xn_zone(ih1)  = 0.75_rt - HALF*metalicity
     xn_zone(ihe4) = 0.25_rt - HALF*metalicity

     do jj = lo(2), hi(2)
        temp_zone = 10.0_rt**(log10(temp_min) + dble(jj)*dlogT)

        do ii = lo(1), hi(1)
           dens_zone = 10.0_rt**(log10(dens_min) + dble(ii)*dlogrho)

           eos_state % rho = dens_zone
           eos_state % T = temp_zone
           eos_state % xn(:) = xn_zone(:)

           ! store default state
           sp(ii, jj, kk, p % irho) = dens_zone
           sp(ii, jj, kk, p % itemp) = temp_zone
           sp(ii, jj, kk, p % ispec: p % ispec-1+nspec) = xn_zone(:)

           ! call EOS using rho, T
           call eos(eos_input_rt, eos_state)

           ! call the conductivity routine
           call actual_conductivity(eos_state)

           sp(ii, jj, kk, p % ih) = eos_state % h
           sp(ii, jj, kk, p % ie) = eos_state % e
           sp(ii, jj, kk, p % ip) = eos_state % p
           sp(ii, jj, kk, p % is) = eos_state % s
           
           sp(ii, jj, kk, p % iconductivity) = eos_state % conductivity
        enddo
     enddo
  enddo

end subroutine do_conductivity
