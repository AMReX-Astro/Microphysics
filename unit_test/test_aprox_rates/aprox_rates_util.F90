subroutine do_rates(lo, hi, &
                    sp, s_lo, s_hi, npts) bind(C, name="do_rates")

    use variables
    use network
    use eos_type_module
    use eos_module
    use amrex_fort_module, only : rt => amrex_real
    use amrex_constants_module
    use extern_probin_module
    use aprox_rates_module

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: sp(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), p % n_plot_comps)
    integer, intent(in), value :: npts

    real(rt) :: dlogrho, dlogT, dNi
    real(rt) :: ni56, temp_zone, dens_zone
    real(rt) :: xn_zone(nspec)
    real(rt) :: fr, dfrdt, rr, drrdt, dfrdd, drrdd
    real(rt) :: rn56ec, sn56ec ! langanke
    real(rt) :: rpen, rnep, spenc, snepc ! ecapnuc

    type(eos_t) :: eos_state
    type(tf_t) :: tf

    integer :: i, j, k

    !$gpu

    dlogrho = (log10(dens_max) - log10(dens_min))/(npts - 1)
    dlogT   = (log10(temp_max) - log10(temp_min))/(npts - 1)
    dNi    = ONE/(npts - 1)

    do k = lo(3), hi(3)

        ! set the composition -- approximately solar
        ni56 = ZERO + dble(k)*dNi
        xn_zone(:) = ONE - ni56/(nspec - 1)   ! all but H, He
        xn_zone(ini56)  = ni56

        do j = lo(2), hi(2)
            temp_zone = 10.0_rt**(log10(temp_min) + dble(j)*dlogT)

            call get_tfactors(temp_zone, tf)

            do i = lo(1), hi(1)
                dens_zone = 10.0_rt**(log10(dens_min) + dble(i)*dlogrho)

                eos_state % rho = dens_zone
                eos_state % T = temp_zone
                eos_state % xn(:) = xn_zone(:)

                ! store default state
                sp(i, j, k, p % irho) = dens_zone
                sp(i, j, k, p % itemp) = temp_zone
                sp(i, j, k, p % ini56) = ni56

                ! call EOS using rho, T
                call eos(eos_input_rt, eos_state)

                call rate_c12ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ic12ag) = fr
                sp(i, j, k, p % ic12ag+1) = dfrdt
                sp(i, j, k, p % ic12ag+2) = rr
                sp(i, j, k, p % ic12ag+3) = drrdt

                call rate_c12ag_deboer17(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ic12ag_deboer17) = fr
                sp(i, j, k, p % ic12ag_deboer17+1) = dfrdt
                sp(i, j, k, p % ic12ag_deboer17+2) = rr
                sp(i, j, k, p % ic12ag_deboer17+3) = drrdt

                call rate_tripalf(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % itriplealf) = fr
                sp(i, j, k, p % itriplealf+1) = dfrdt
                sp(i, j, k, p % itriplealf+2) = rr
                sp(i, j, k, p % itriplealf+3) = drrdt

                call rate_c12c12(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ic12c12) = fr
                sp(i, j, k, p % ic12c12+1) = dfrdt
                sp(i, j, k, p % ic12c12+2) = rr
                sp(i, j, k, p % ic12c12+3) = drrdt

                call rate_c12o16(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ic12o16) = fr
                sp(i, j, k, p % ic12o16+1) = dfrdt
                sp(i, j, k, p % ic12o16+2) = rr
                sp(i, j, k, p % ic12o16+3) = drrdt

                call rate_o16o16(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % io16o16) = fr
                sp(i, j, k, p % io16o16+1) = dfrdt
                sp(i, j, k, p % io16o16+2) = rr
                sp(i, j, k, p % io16o16+3) = drrdt

                call rate_o16ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % io16ag) = fr
                sp(i, j, k, p % io16ag+1) = dfrdt
                sp(i, j, k, p % io16ag+2) = rr
                sp(i, j, k, p % io16ag+3) = drrdt

                call rate_ne20ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ine20ag) = fr
                sp(i, j, k, p % ine20ag+1) = dfrdt
                sp(i, j, k, p % ine20ag+2) = rr
                sp(i, j, k, p % ine20ag+3) = drrdt

                call rate_mg24ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % img24ag) = fr
                sp(i, j, k, p % img24ag+1) = dfrdt
                sp(i, j, k, p % img24ag+2) = rr
                sp(i, j, k, p % img24ag+3) = drrdt

                call rate_mg24ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % img24ap) = fr
                sp(i, j, k, p % img24ap+1) = dfrdt
                sp(i, j, k, p % img24ap+2) = rr
                sp(i, j, k, p % img24ap+3) = drrdt

                call rate_al27pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ial27pg) = fr
                sp(i, j, k, p % ial27pg+1) = dfrdt
                sp(i, j, k, p % ial27pg+2) = rr
                sp(i, j, k, p % ial27pg+3) = drrdt

                call rate_al27pg_old(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ial27pg_old) = fr
                sp(i, j, k, p % ial27pg_old+1) = dfrdt
                sp(i, j, k, p % ial27pg_old+2) = rr
                sp(i, j, k, p % ial27pg_old+3) = drrdt

                call rate_si28ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % isi28ag) = fr
                sp(i, j, k, p % isi28ag+1) = dfrdt
                sp(i, j, k, p % isi28ag+2) = rr
                sp(i, j, k, p % isi28ag+3) = drrdt
                
                call rate_si28ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % isi28ap) = fr
                sp(i, j, k, p % isi28ap+1) = dfrdt
                sp(i, j, k, p % isi28ap+2) = rr
                sp(i, j, k, p % isi28ap+3) = drrdt

                call rate_p31pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ip31pg) = fr
                sp(i, j, k, p % ip31pg+1) = dfrdt
                sp(i, j, k, p % ip31pg+2) = rr
                sp(i, j, k, p % ip31pg+3) = drrdt

                call rate_s32ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % is32ag) = fr
                sp(i, j, k, p % is32ag+1) = dfrdt
                sp(i, j, k, p % is32ag+2) = rr
                sp(i, j, k, p % is32ag+3) = drrdt

                call rate_s32ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % is32ap) = fr
                sp(i, j, k, p % is32ap+1) = dfrdt
                sp(i, j, k, p % is32ap+2) = rr
                sp(i, j, k, p % is32ap+3) = drrdt

                call rate_cl35pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % icl35pg) = fr
                sp(i, j, k, p % icl35pg+1) = dfrdt
                sp(i, j, k, p % icl35pg+2) = rr
                sp(i, j, k, p % icl35pg+3) = drrdt

                call rate_ar36ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % iar36ag) = fr
                sp(i, j, k, p % iar36ag+1) = dfrdt
                sp(i, j, k, p % iar36ag+2) = rr
                sp(i, j, k, p % iar36ag+3) = drrdt

                call rate_ar36ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % iar36ap) = fr
                sp(i, j, k, p % iar36ap+1) = dfrdt
                sp(i, j, k, p % iar36ap+2) = rr
                sp(i, j, k, p % iar36ap+3) = drrdt

                call rate_k39pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ik39pg) = fr
                sp(i, j, k, p % ik39pg+1) = dfrdt
                sp(i, j, k, p % ik39pg+2) = rr
                sp(i, j, k, p % ik39pg+3) = drrdt

                call rate_ca40ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ica40ag) = fr
                sp(i, j, k, p % ica40ag+1) = dfrdt
                sp(i, j, k, p % ica40ag+2) = rr
                sp(i, j, k, p % ica40ag+3) = drrdt

                call rate_ca40ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ica40ap) = fr
                sp(i, j, k, p % ica40ap+1) = dfrdt
                sp(i, j, k, p % ica40ap+2) = rr
                sp(i, j, k, p % ica40ap+3) = drrdt

                call rate_sc43pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % isc43pg) = fr
                sp(i, j, k, p % isc43pg+1) = dfrdt
                sp(i, j, k, p % isc43pg+2) = rr
                sp(i, j, k, p % isc43pg+3) = drrdt

                call rate_ti44ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % iti44ag) = fr
                sp(i, j, k, p % iti44ag+1) = dfrdt
                sp(i, j, k, p % iti44ag+2) = rr
                sp(i, j, k, p % iti44ag+3) = drrdt

                call rate_ti44ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % iti44ap) = fr
                sp(i, j, k, p % iti44ap+1) = dfrdt
                sp(i, j, k, p % iti44ap+2) = rr
                sp(i, j, k, p % iti44ap+3) = drrdt

                call rate_v47pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % iv47pg) = fr
                sp(i, j, k, p % iv47pg+1) = dfrdt
                sp(i, j, k, p % iv47pg+2) = rr
                sp(i, j, k, p % iv47pg+3) = drrdt

                call rate_cr48ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % icr48ag) = fr
                sp(i, j, k, p % icr48ag+1) = dfrdt
                sp(i, j, k, p % icr48ag+2) = rr
                sp(i, j, k, p % icr48ag+3) = drrdt

                call rate_cr48ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % icr48ap) = fr
                sp(i, j, k, p % icr48ap+1) = dfrdt
                sp(i, j, k, p % icr48ap+2) = rr
                sp(i, j, k, p % icr48ap+3) = drrdt

                call rate_mn51pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % imn51pg) = fr
                sp(i, j, k, p % imn51pg+1) = dfrdt
                sp(i, j, k, p % imn51pg+2) = rr
                sp(i, j, k, p % imn51pg+3) = drrdt

                call rate_fe52ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife52ag) = fr
                sp(i, j, k, p % ife52ag+1) = dfrdt
                sp(i, j, k, p % ife52ag+2) = rr
                sp(i, j, k, p % ife52ag+3) = drrdt

                call rate_fe52ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife52ap) = fr
                sp(i, j, k, p % ife52ap+1) = dfrdt
                sp(i, j, k, p % ife52ap+2) = rr
                sp(i, j, k, p % ife52ap+3) = drrdt

                call rate_co55pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ico55pg) = fr
                sp(i, j, k, p % ico55pg+1) = dfrdt
                sp(i, j, k, p % ico55pg+2) = rr
                sp(i, j, k, p % ico55pg+3) = drrdt

                call rate_pp(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ipp) = fr
                sp(i, j, k, p % ipp+1) = dfrdt
                sp(i, j, k, p % ipp+2) = rr
                sp(i, j, k, p % ipp+3) = drrdt

                call rate_png(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ipng) = fr
                sp(i, j, k, p % ipng+1) = dfrdt
                sp(i, j, k, p % ipng+2) = rr
                sp(i, j, k, p % ipng+3) = drrdt

                call rate_dpg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % idpg) = fr
                sp(i, j, k, p % idpg+1) = dfrdt
                sp(i, j, k, p % idpg+2) = rr
                sp(i, j, k, p % idpg+3) = drrdt

                call rate_he3ng(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ihe3ng) = fr
                sp(i, j, k, p % ihe3ng+1) = dfrdt
                sp(i, j, k, p % ihe3ng+2) = rr
                sp(i, j, k, p % ihe3ng+3) = drrdt

                call rate_he3he3(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ihe3he3) = fr
                sp(i, j, k, p % ihe3he3+1) = dfrdt
                sp(i, j, k, p % ihe3he3+2) = rr
                sp(i, j, k, p % ihe3he3+3) = drrdt

                call rate_he3he4(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ihe3he4) = fr
                sp(i, j, k, p % ihe3he4+1) = dfrdt
                sp(i, j, k, p % ihe3he4+2) = rr
                sp(i, j, k, p % ihe3he4+3) = drrdt

                call rate_c12pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ic12pg) = fr
                sp(i, j, k, p % ic12pg+1) = dfrdt
                sp(i, j, k, p % ic12pg+2) = rr
                sp(i, j, k, p % ic12pg+3) = drrdt

                call rate_n14pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % in14pg) = fr
                sp(i, j, k, p % in14pg+1) = dfrdt
                sp(i, j, k, p % in14pg+2) = rr
                sp(i, j, k, p % in14pg+3) = drrdt

                call rate_n15pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % in15pg) = fr
                sp(i, j, k, p % in15pg+1) = dfrdt
                sp(i, j, k, p % in15pg+2) = rr
                sp(i, j, k, p % in15pg+3) = drrdt

                call rate_n15pa(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % in15pa) = fr
                sp(i, j, k, p % in15pa+1) = dfrdt
                sp(i, j, k, p % in15pa+2) = rr
                sp(i, j, k, p % in15pa+3) = drrdt

                call rate_o16pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % io16pg) = fr
                sp(i, j, k, p % io16pg+1) = dfrdt
                sp(i, j, k, p % io16pg+2) = rr
                sp(i, j, k, p % io16pg+3) = drrdt

                call rate_n14ag(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % in14ag) = fr
                sp(i, j, k, p % in14ag+1) = dfrdt
                sp(i, j, k, p % in14ag+2) = rr
                sp(i, j, k, p % in14ag+3) = drrdt

                call rate_fe52ng(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife52ng) = fr
                sp(i, j, k, p % ife52ng+1) = dfrdt
                sp(i, j, k, p % ife52ng+2) = rr
                sp(i, j, k, p % ife52ng+3) = drrdt

                call rate_fe53ng(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife53ng) = fr
                sp(i, j, k, p % ife53ng+1) = dfrdt
                sp(i, j, k, p % ife53ng+2) = rr
                sp(i, j, k, p % ife53ng+3) = drrdt

                call rate_fe54ng(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife54ng) = fr
                sp(i, j, k, p % ife54ng+1) = dfrdt
                sp(i, j, k, p % ife54ng+2) = rr
                sp(i, j, k, p % ife54ng+3) = drrdt

                call rate_fe54pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife54pg) = fr
                sp(i, j, k, p % ife54pg+1) = dfrdt
                sp(i, j, k, p % ife54pg+2) = rr
                sp(i, j, k, p % ife54pg+3) = drrdt

                call rate_fe54ap(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife54ap) = fr
                sp(i, j, k, p % ife54ap+1) = dfrdt
                sp(i, j, k, p % ife54ap+2) = rr
                sp(i, j, k, p % ife54ap+3) = drrdt

                call rate_fe55ng(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife55ng) = fr
                sp(i, j, k, p % ife55ng+1) = dfrdt
                sp(i, j, k, p % ife55ng+2) = rr
                sp(i, j, k, p % ife55ng+3) = drrdt

                call rate_fe56pg(tf, dens_zone, fr, dfrdt, dfrdd, rr, drrdt, drrdd)

                sp(i, j, k, p % ife56pg) = fr
                sp(i, j, k, p % ife56pg+1) = dfrdt
                sp(i, j, k, p % ife56pg+2) = rr
                sp(i, j, k, p % ife56pg+3) = drrdt

                call langanke(temp_zone, dens_zone, eos_state % xn(ini56), &
                        eos_state % y_e, rn56ec, sn56ec)

                sp(i, j, k, p % ilanganke) = rn56ec
                sp(i, j, k, p % ilanganke+1) = sn56ec

                call ecapnuc(eos_state % eta, dens_zone, rpen, rnep, spenc, snepc)

                sp(i, j, k, p % iecapnuc) = rpen
                sp(i, j, k, p % iecapnuc+1) = rnep
                sp(i, j, k, p % iecapnuc+2) = spenc
                sp(i, j, k, p % iecapnuc+3) = snepc

            enddo
        enddo
    enddo

end subroutine do_rates
