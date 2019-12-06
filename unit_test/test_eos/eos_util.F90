subroutine do_eos(lo, hi, &
                  sp, s_lo, s_hi, npts) bind(C, name="do_eos")

  use variables
  use network
  use eos_type_module, only: eos_t, eos_input_rt, eos_input_rp, eos_input_tp, eos_input_re, &
                             eos_input_ph, eos_input_ps, eos_input_rh, eos_input_th, &
                             copy_eos_t
  use eos_module, only : eos
  use eos_composition_module, only : eos_xderivs_t, composition, composition_derivatives
  use microphysics_type_module
  use extern_probin_module
  use actual_eos_module, only : eos_name

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
  type(eos_xderivs_t) :: eos_xderivs

  integer :: ii, jj, kk, n

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
           call composition_derivatives(eos_state, eos_xderivs)

           call copy_eos_t(eos_state_reference, eos_state)

           sp(ii, jj, kk, p % ih) = eos_state % h
           sp(ii, jj, kk, p % ie) = eos_state % e
           sp(ii, jj, kk, p % ip) = eos_state % p
           sp(ii, jj, kk, p % is) = eos_state % s

           sp(ii, jj, kk, p % icv) = eos_state % cv
           sp(ii, jj, kk, p % icp) = eos_state % cp
           sp(ii, jj, kk, p % ine) = eos_state % xne
           sp(ii, jj, kk, p % inp) = eos_state % xnp
           sp(ii, jj, kk, p % ieta) = eos_state % eta
           sp(ii, jj, kk, p % ipele) = eos_state % pele
           sp(ii, jj, kk, p % ippos) = eos_state % ppos
           sp(ii, jj, kk, p % imu) = eos_state % mu
           sp(ii, jj, kk, p % imue) = eos_state % mu_e
           sp(ii, jj, kk, p % idpdt) = eos_state % dpdt
           sp(ii, jj, kk, p % idpdr) = eos_state % dpdr
           sp(ii, jj, kk, p % idedt) = eos_state % dedt
           sp(ii, jj, kk, p % idedr) = eos_state % dedr
           sp(ii, jj, kk, p % idhdt) = eos_state % dhdt
           sp(ii, jj, kk, p % idhdr) = eos_state % dhdr
           sp(ii, jj, kk, p % idsdt) = eos_state % dsdt
           sp(ii, jj, kk, p % idsdr) = eos_state % dsdr
           do n = 0, nspec-1
              sp(ii, jj, kk, p % idpdx + n) = eos_xderivs % dpdx(n+1)
              sp(ii, jj, kk, p % idedx + n) = eos_xderivs % dedx(n+1)
              sp(ii, jj, kk, p % idhdx + n) = eos_xderivs % dhdx(n+1)
           end do
           sp(ii, jj, kk, p % igam1) = eos_state % gam1
           sp(ii, jj, kk, p % ics) = eos_state % cs
           sp(ii, jj, kk, p % iabar) = eos_state % abar
           sp(ii, jj, kk, p % izbar) = eos_state % zbar
           sp(ii, jj, kk, p % idpda) = eos_state % dpda
           sp(ii, jj, kk, p % idpdz) = eos_state % dpdz
           sp(ii, jj, kk, p % ideda) = eos_state % deda
           sp(ii, jj, kk, p % idedz) = eos_state % dedz
           sp(ii, jj, kk, p % idpde) = eos_state % dpde
           sp(ii, jj, kk, p % idpdre) = eos_state % dpdr_e


           ! call EOS using rho, h

           ! reset T to give it some work to do
           eos_state % T = 100.e0_rt

           call eos(eos_input_rh, eos_state)

           sp(ii, jj, kk, p % ierr_T_eos_rh) = &
                abs(eos_state % T - temp_zone)/temp_zone

           call copy_eos_t(eos_state, eos_state_reference)

           
           ! call EOS using T, p

           ! reset rho to give it some work to do
           eos_state % rho = 1.e0_rt

           call eos(eos_input_tp, eos_state)

           sp(ii, jj, kk, p % ierr_rho_eos_tp) = &
                abs(eos_state % rho - dens_zone)/dens_zone

           call copy_eos_t(eos_state, eos_state_reference)

           
           ! call EOS using r, p

           ! reset T to give it some work to do
           eos_state % T = 100.e0_rt

           call eos(eos_input_rp, eos_state)

           sp(ii, jj, kk, p % ierr_T_eos_rp) = &
                abs(eos_state % T - temp_zone)/temp_zone

           call copy_eos_t(eos_state, eos_state_reference)


           ! call EOS using r, e

           ! reset T to give it some work to do
           eos_state % T = 100.e0_rt

           call eos(eos_input_re, eos_state)

           sp(ii, jj, kk, p % ierr_T_eos_re) = &
                abs(eos_state % T - temp_zone)/temp_zone

           call copy_eos_t(eos_state, eos_state_reference)


           ! call EOS using p, s

           ! reset T and rho to give it some work to do
           eos_state % T = 100.e0_rt
           eos_state % rho = 1.e0_rt


           ! some EOSes don't have physically valid treatments
           ! of entropy throughout the entire rho-T plane
           if (eos_state%s > ZERO) then

              call eos(eos_input_ps, eos_state)

              ! store the thermodynamic state
              sp(ii, jj, kk, p % ierr_T_eos_ps) = &
                   abs(eos_state % T - temp_zone)/temp_zone
              sp(ii, jj, kk, p % ierr_rho_eos_ps) = &
                   abs(eos_state % rho - dens_zone)/dens_zone

           else
              sp(ii, jj, kk, p % ierr_T_eos_ps) = ZERO
              sp(ii, jj, kk, p % ierr_rho_eos_ps) = ZERO

           endif

           call copy_eos_t(eos_state, eos_state_reference)


           ! call EOS using p, h

           ! reset T and rho to give it some work to do
           eos_state % T = 100.e0_rt
           eos_state % rho = 1.e0_rt

           call eos(eos_input_ph, eos_state)

           sp(ii, jj, kk, p % ierr_T_eos_ph) = &
                abs(eos_state % T - temp_zone)/temp_zone
           sp(ii, jj, kk, p % ierr_rho_eos_ph) = &
                abs(eos_state % rho - dens_zone)/dens_zone

           call copy_eos_t(eos_state, eos_state_reference)


           ! call EOS using T, h
           ! this doesn't work for all EOSes (where h doesn't depend on T)
#ifndef EOS_GAMMA_LAW_GENERAL
           ! reset rho to give it some work to do -- for helmeos, h is not
           ! monotonic, so we only perturb rho slightly here
           eos_state % rho = 0.9_rt * eos_state % rho

           call eos(eos_input_th, eos_state)

           sp(ii, jj, kk, p % ierr_rho_eos_th) = &
                abs(eos_state % rho - dens_zone)/dens_zone

           call copy_eos_t(eos_state, eos_state_reference)

#else
           sp(ii, jj, kk, p % ierr_rho_eos_th) = ZERO
#endif
        enddo
     enddo
  enddo

end subroutine do_eos
