!***************************************************************************************************
! output.f90 10/18/17
! This file contains the ts_output routine, which controls the information output at the end of each
! timestep, and the final_output routine, which controls the output at the end of the evolution.
!***************************************************************************************************

Subroutine ts_output(kstep,enuc,edot,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! The per timestep output routine. If the flag itsout > 0, full_net calls this routine to handle
  ! stepwise output.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: y
  Use conditions, Only: rho, t, t9, tdel
  Use controls, Only: nnucout, idiag, inucout, itsout, kmon, lun_diag, lun_ev, lun_stdout, lun_ts, &
    & nzbatchmx, szbatch, lzactive
  Use flux_data, Only: flx, flx_int
  Use match_data, Only: iwflx, mflx, nflx
  Use nuc_number, Only: ny
  Use nuclear_data, Only: aa, nname
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_output
  Use xnet_interface, Only: flux
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep
  Real(dp), Intent(in) :: enuc(:), edot(:)

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Character(40) :: ev_format
  Integer :: i, j, k, izb, izone
  Logical, Pointer :: mask(:)

  start_timer = xnet_wtime()
  timer_output = timer_output - start_timer

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  ! Calculate reaction fluxes
  If ( itsout > 0 .or. idiag > 0 ) Then
    If ( kstep > 0 ) Then
      Call flux(mask_in = mask)
    Else
      flx_int = 0.0
    EndIf
  EndIf

  If ( itsout > 0 ) Then
    Write(ev_format,"(a28,i2,a10)") "(i4,1es15.8,2es10.3,2es10.2,",nnucout,"es9.2,2i2)"
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1

        ! An abundance snapshot is written to the binary file
        Write(lun_ts(izb)) kstep,t(izb),t9(izb),rho(izb),tdel(izb),edot(izb),y(:,izb),flx(:,izb)

        ! For itsout>=2, output important mass fractions to the ASCII file
        If( itsout >= 2 ) Write(lun_ev(izb),ev_format) &
          & kstep,t(izb),t9(izb),rho(izb),edot(izb),tdel(izb), &
          & (aa(inucout(i))*y(inucout(i),izb),i=1,nnucout),(kmon(j,izb),j=1,2)

        ! For itsout>=3, output time and thermo evolution to the screen
        If( itsout >= 3 ) Write(lun_stdout,"(3i5,4es12.4,2i3)") &
          & izb, izone, kstep,t(izb),tdel(izb),t9(izb),rho(izb),(kmon(j,izb),j=1,2)

        ! For itsout>=4, output abundances for each timestep to the diagnostic file
        If( itsout >= 4 ) Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i,izb),i=1,ny)

        ! For itsout>=5, output fluxes for each timestep to the diagnostic file
        If( itsout >= 5 ) Write(lun_diag,"(i5,8a5,i5,es11.3)") &
          & (k,(nname(nflx(j,k)),j=1,3),' <-> ',(nname(nflx(j,k)),j=4,7),iwflx(k),flx(k,izb),k=1,mflx)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_output = timer_output + stop_timer

  Return
End Subroutine ts_output

Subroutine final_output(kstep,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! full_net calls this routine to handle output at the end of its execution.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: y
  Use conditions, Only: tt
  Use controls, Only: changemx, iconvc, isolv, kitmx, kstmx, ktot, lun_diag, tolc, tolm, yacc, &
    & nzbatchmx, szbatch, lzactive, idiag
  Use flux_data, Only: flx_int
  Use match_data, Only: iwflx, mflx, nflx
  Use nuc_number, Only: ny
  Use nuclear_data, Only: nname, aa
  Use thermo_data, Only: tstop
  Use timers
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer, Parameter :: itimer_reset = 0
  Integer :: iflx_mx
  Integer :: i, k, izb, izone
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  If ( idiag > 0 ) Then

    start_timer = xnet_wtime()
    timer_output = timer_output - start_timer

    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1

        ! Write final abundances to diagnostic output (in ASCII)
        Write(lun_diag,"(a3,3i6,2es14.7)") 'End',izone,kstep,kstmx,tt(izb),tstop(izb)
        Write(lun_diag,"(4(a5,es14.7,1x))") (nname(i),aa(i)*y(i,izb),i=1,ny)
        Write(lun_diag,"(a3,3i6,4es10.3)") 'CNT',isolv,iconvc,kitmx,yacc,tolm,tolc,changemx

        ! Write integrated fluxes to diagnotic output
        iflx_mx = maxloc(abs(flx_int(:,izb)),dim=1)
        Write(lun_diag,"(2x,a8,i5,es11.3,8a5)") 'Flux Max',&
          & iflx_mx,flx_int(iflx_mx,izb),(nname(nflx(i,iflx_mx)),i=1,3),' <-> ',(nname(nflx(i,iflx_mx)),i=4,7)
        Write(lun_diag,"(i5,8a5,i5,es11.3)") &
          & (k,(nname(nflx(i,k)),i=1,3),' <-> ',(nname(nflx(i,k)),i=4,7),iwflx(k),flx_int(k,izb),k=1,mflx)

        Write(lun_diag,"(a10,a5,2a10)") ' Counters: ','Zone','Timestep','Iteration'
        Write(lun_diag,"(10x,i5,2i10)") izone,ktot(1,izb),ktot(2,izb)
      EndIf
    EndDo

    stop_timer = xnet_wtime()
    timer_output = timer_output + stop_timer

    ! Write timers
    Write(lun_diag,"(a8,8a10)") 'Timers: ','Total','Solver','Jacobian','Deriv','CrossSect','Screening','Setup','Output'
    Write(lun_diag,"(8x,8es10.3)") timer_xnet,timer_solve,timer_jacob,timer_deriv,timer_csect,timer_scrn,timer_setup,timer_output

    ! Use the following line to restart the timers
    If ( itimer_reset > 0 ) Call reset_timers

  End If

  Return
End Subroutine final_output