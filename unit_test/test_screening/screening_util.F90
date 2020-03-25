subroutine do_screening(lo, hi, &
                        sp, s_lo, s_hi, npts) bind(C, name="do_screening")

  use variables
  use network
  use eos_type_module
  use eos_module
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use extern_probin_module
  use screening_module, only: add_screening_factor
  use screening_module, only: screen5, plasma_state, fill_plasma_state

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: s_lo(3), s_hi(3)
  real(rt), intent(inout) :: sp(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), p % n_plot_comps)
  integer, intent(in), value :: npts

  real(rt) :: dlogrho, dlogT, dmetal
  real(rt) :: metalicity, temp_zone, dens_zone
  real(rt) :: xn_zone(nspec), ymass(nspec)

  real(rt) :: sc1a,sc1adt,sc1add

  type (plasma_state) :: state

  integer :: ii, jj, kk

  integer :: jscr

  !$gpu

  ! initialize the screening factors
  call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ihe4),aion(ihe4),4.0e0_rt,8.0e0_rt)

  call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))

  call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))

  call add_screening_factor(zion(io16),aion(io16),zion(io16),aion(io16))

  call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))

  call add_screening_factor(13.0e0_rt,27.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(isi28),aion(isi28),zion(ihe4),aion(ihe4))

  call add_screening_factor(15.0e0_rt,31.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(is32),aion(is32),zion(ihe4),aion(ihe4))

  call add_screening_factor(17.0e0_rt,35.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(iar36),aion(iar36),zion(ihe4),aion(ihe4))

  call add_screening_factor(19.0e0_rt,39.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ica40),aion(ica40),zion(ihe4),aion(ihe4))

  call add_screening_factor(21.0e0_rt,43.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(iti44),aion(iti44),zion(ihe4),aion(ihe4))

  call add_screening_factor(23.0e0_rt,47.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(icr48),aion(icr48),zion(ihe4),aion(ihe4))

  call add_screening_factor(25.0e0_rt,51.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ife52),aion(ife52),zion(ihe4),aion(ihe4))

  call add_screening_factor(27.0e0_rt,55.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ife54),aion(ife54),1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ife54),aion(ife54),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ife56),aion(ife56),1.0e0_rt,1.0e0_rt)

  call add_screening_factor(1.0e0_rt,2.0e0_rt,zion(ih1),aion(ih1))

  call add_screening_factor(zion(ih1),aion(ih1),zion(ih1),aion(ih1))

  call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3))

  call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ic12),aion(ic12),zion(ih1),aion(ih1))

  call add_screening_factor(zion(in14),aion(in14),zion(ih1),aion(ih1))

  call add_screening_factor(zion(io16),aion(io16),zion(ih1),aion(ih1))

  call add_screening_factor(zion(in14),aion(in14),zion(ihe4),aion(ihe4))


  dlogrho = (log10(dens_max) - log10(dens_min))/(npts - 1)
  dlogT   = (log10(temp_max) - log10(temp_min))/(npts - 1)
  dmetal    = (metalicity_max  - ZERO)/(npts - 1)

  do kk = lo(3), hi(3)

     ! set the composition -- approximately solar
     metalicity = ZERO + dble(kk)*dmetal
     xn_zone(:) = metalicity/(nspec - 2)   ! all but H, He
     xn_zone(ih1)  = 0.75_rt - HALF*metalicity
     xn_zone(ihe4) = 0.25_rt - HALF*metalicity

     ymass(:) = xn_zone(:) / aion(:)

     do jj = lo(2), hi(2)
        temp_zone = 10.0_rt**(log10(temp_min) + dble(jj)*dlogT)

        do ii = lo(1), hi(1)
           dens_zone = 10.0_rt**(log10(dens_min) + dble(ii)*dlogrho)

           ! store default state
           sp(ii, jj, kk, p % irho) = dens_zone
           sp(ii, jj, kk, p % itemp) = temp_zone
           sp(ii, jj, kk, p % ispec: p % ispec-1+nspec) = xn_zone(:)

           call fill_plasma_state(state, temp_zone, dens_zone, ymass(1:nspec))


           ! 3-alpha
           jscr = 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_he4_he4) = sc1a

           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_he4_be8) = sc1a

           ! c12(a,g)o16
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_c12_he4) = sc1a

           ! c12 + c12
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_c12_c12) = sc1a

           ! c12 + o16
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_c12_o16) = sc1a

           ! o16 + o16
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_o16_o16) = sc1a

           ! o16 + ne20
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_o16_he4) = sc1a

           ! ne20(a,g)mg24
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_ne20_he4) = sc1a

           ! mg24(a,g)si28
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_mg24_he4) = sc1a

           ! al27(p,g)si28
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_al27_p) = sc1a

           ! si28 tp s32
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_si28_he4) = sc1a

           ! p31(p,g)s32
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_p31_p) = sc1a

           ! s32 to ar36
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_s32_he4) = sc1a

           ! cl35(p,g)ar36
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_cl35_p) = sc1a

           ! ar36 to ca40
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_ar36_he4) = sc1a

           ! k39(p,g)ca40
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_k39_p) = sc1a

           ! ca40 to ti44
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_ca40_he4) = sc1a

           ! sc43(p,g)ti44
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_sc43_p) = sc1a

           ! ti44 to cr48
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_ti44_he4) = sc1a

           ! v47(p,g)cr48
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_v47_p) = sc1a

           ! cr48 to fe52
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_cr48_he4) = sc1a

           ! mn51(p,g)fe52
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_mn51_p) = sc1a

           ! fe to ni
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_fe52_he4) = sc1a

           ! co55(p,g)ni56
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_co55_p) = sc1a

           ! fe54(p,g)co55
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_fe54_p) = sc1a

           ! fe54(a,p)co57
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_fe54_he4) = sc1a

           ! fe56(p,g)co57
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_fe56_p) = sc1a

           ! d + p
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_d_p) = sc1a

           ! pp
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_p_p) = sc1a

           ! he3 + he3
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_he3_he3) = sc1a

           ! he3 + he4
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_he3_he4) = sc1a

           ! c12(p,g)n13
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_c12_p) = sc1a

           ! n14(p,g)o15
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_n14_p) = sc1a

           ! o16(p,g)f17
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_o16_p) = sc1a

           ! n14(a,g)f18
           jscr = jscr + 1
           call screen5(state,jscr,sc1a,sc1adt,sc1add)
           sp(ii, jj, kk, p % iscn_n14_he4) = sc1a

        enddo
     enddo
  enddo

end subroutine do_screening
