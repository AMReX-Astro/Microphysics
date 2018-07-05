subroutine do_eos(lo, hi, &
                  state, s_lo, s_hi)

  use variables

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), p % n_plot_comps)

  nrho = extent(mla%mba%pd(1),1)
  nT = extent(mla%mba%pd(1),2)
  nX = extent(mla%mba%pd(1),3)

  dlogrho = (log10(dens_max) - log10(dens_min))/(nrho - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(nT - 1)
  dmetal    = (metalicity_max  - ZERO)/(nX - 1)

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


end program test_react
