!***************************************************************************************************
! full_net.f90 10/18/17
! The subroutines in this file perform nucleosynthesis for a single thermodynamic trajectory
! (usually a Lagrangian mass particle or zone). The handling of multiple trajectories, including
! multi-processing, is carried out externally.
!***************************************************************************************************

Subroutine full_net(kstep)
  !-------------------------------------------------------------------------------------------------
  ! The abundance evolution is performed over a series of timesteps, with the duration of the
  ! timestep determined by the integration scheme and the changing thermodynamic conditions.
  ! Integration is performed by a choice of methods controlled by the isolv flag.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: yo, y, yt, ydot
  Use conditions, Only: t, to, tt, tdel, tdel_old, t9, t9o, t9t, t9dot, rho, rhoo, rhot, yeo, ye, yet, &
    nt, nto, ntt
  Use controls
  Use nuc_number, Only: ny
  Use nuclear_data, Only: nname
  Use thermo_data, Only: tstart, tstop
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_xnet
  Use xnet_interface, Only: benuc, final_output, solve_be, t9rhofind, timestep, ts_output, &
    update_iweak, xnet_terminate
  Use xnet_types, Only: dp
  Implicit None

  ! Output variables
  Integer, Intent(out) :: kstep

  ! Local variables
  Integer, Parameter :: kstep_output = 10
  Real(dp) :: enm(nzbatchmx), enb(nzbatchmx), enold(nzbatchmx), en0(nzbatchmx), edot(nzbatchmx)
  Real(dp) :: ytot, ztot, atot
  Integer :: idiag0, its(nzbatchmx), mykstep(nzbatchmx)
  Integer :: i, j, k, izb, izone
  Logical :: lzstep(nzbatchmx)

  start_timer = xnet_wtime()
  timer_xnet = timer_xnet - start_timer

  ! Initialize counters
  kstep = 0
  kmon = 0
  ktot = 0

  ! Set reaction controls not read in from control
  iweak = iweak0
  idiag0 = idiag

  ! Initialize trial time step abundances and conditions
  Call t9rhofind(0,t,nt,t9,rho)
  Call update_iweak(t9)
  nto = nt   ; ntt = nt
  to = t     ; tt = t
  yo = y     ; yt = y
  t9o = t9   ; t9t = t9
  rhoo = rho ; rhot = rho
  yet = ye   ; yeo = ye
  tdel_old = tdel

  edot = 0.0 ; en0 = 0.0 ; enm = 0.0
  Do izb = 1, nzbatchmx
    If ( lzactive(izb) ) Then

      ! Calculate the total energy of the nuclei
      Call benuc(yt(:,izb),enb(izb),enm(izb),ytot,ztot,atot)
      en0(izb) = enm(izb)
      edot(izb) = 0.0

      If ( itsout > 0 ) Write(lun_stdout,"(a,i6,a,i2,2(a,es10.3))") &
        & 'Max Step',kstmx,' IDiag=',idiag,' Start Time',tstart(izb),' Stop Time',tstop(izb)
    EndIf
  EndDo

  ! Output initial abundances and conditions
  Call ts_output(kstep,enm-en0,edot)

  ! Start evolution
  Where ( lzactive )
    its = 0
  ElseWhere
    its = -1
  EndWhere
  mykstep = 0
  lzstep = ( its < 0 )
  Do kstep = 1, kstmx

    ! Determine if this is an output step
    idiag = idiag0
!   If ( mod(kstep,kstep_output) == 0 ) idiag = 2
    If ( idiag >= 3 ) Write(lun_diag,*) 'KStep',kstep

    ! Calculate an initial guess for the timestep
    Call timestep(kstep,mask_in = (its == 0))
    If ( idiag >= 1 ) Then
      Do izb = 1, nzbatchmx
        If ( its(izb) == 0 ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,i5,2es23.15)") 'TDel',izone,tt(izb),tdel(izb)
        EndIf
      EndDo
    EndIf

    ! Take integration step (only worry about solve_be for now)
!   Select Case (isolv)
!   Case (2)
!     Call solve_bd(kstep,its)
!   Case Default
      Call solve_be(kstep,its)
!   End Select

    Do izb = 1, nzbatchmx

      ! If convergence is successful, output timestep results
      If ( its(izb) == 0 ) Then
        izone = izb + szbatch - 1
        If ( idiag >= 1 ) Then
          Write(lun_diag,"(2i5,5es14.7)") kstep,izone,t(izb),tdel(izb),t9t(izb),rhot(izb),yet(izb)
          If ( idiag >= 2 ) Then
            Write(lun_diag,"(a,es23.15)") 'delta Y',tdel(izb)
            Write(lun_diag,"(a5,4es12.4)") &
              & (nname(k),y(k,izb),yo(k,izb),(y(k,izb)-yo(k,izb)),(tdel(izb)*ydot(k,izb)),k=1,ny)
            If ( iheat > 0 ) Write(lun_diag,"(a,4es12.4)") &
              & 'delta T9',t9(izb),t9o(izb),t9(izb)-t9o(izb),tdel(izb)*t9dot(izb)
          EndIf
!           Write(lun_diag,"(5(a5,es11.4))") (nname(k),y(k,izb),k=1,ny)
        EndIf
        enold(izb) = enm(izb)
        Call benuc(yt(:,izb),enb(izb),enm(izb),ytot,ztot,atot)
        edot(izb) = -(enm(izb)-enold(izb)) / tdel(izb)

      ! If reduced timesteps fail to successfully integrate, warn and flag to remove from loop
      ElseIf ( its(izb) == 1 ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a,i5,a,es12.4,a,es12.4,a)") 'Zone ',izone,' Inter!!!'
        its(izb) = 2
      EndIf
    EndDo
    Call ts_output(kstep,enm-en0,edot,mask_in = (its == 0))

    ! If this zone reaches the stop time, flag it to remove from loop
    Where ( t >= tstop .and. its == 0 )
      mykstep = kstep
      its = -1
    EndWhere

    ! Test if all zones have stopped
    If ( all( its /= 0 ) ) Exit
  EndDo

  ! Test that the stop time is reached
  Do izb = 1, nzbatchmx
    If ( lzactive(izb) ) Then
      If ( t(izb) < tstop(izb) .or. its(izb) > 0 ) Then
        izone = izb + szbatch - 1
        Write(lun_stdout,"(a,i5,a,es12.4,a,es12.4,a)") 'Zone ',izone, &
          & ' Evolution stopped at time=',t(izb),', stop time (',tstop(izb),') not reached!'
        Write(lun_stdout,"(a,i12,a)") &
          & 'Approximately',int((tstop(izb)-t(izb))/tdel(izb)),'more steps needed'
        Call xnet_terminate('[XNet] Evolution failed to converge')
      EndIf
    EndIf
  EndDo

  kstep = max(1,maxval(mykstep))

  stop_timer = xnet_wtime()
  timer_xnet = timer_xnet + stop_timer

  ! Output final state
  Call final_output(kstep)

  Return
End Subroutine full_net
