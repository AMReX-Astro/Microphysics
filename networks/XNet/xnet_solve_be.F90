!***************************************************************************************************
! solve_be.f90 10/18/17
! Backward Euler solver
!
! The routines in this file perform the Backward Euler time integration for the thermonuclear
! reaction network.
!***************************************************************************************************

Subroutine solve_be(kstep,its)
  !-------------------------------------------------------------------------------------------------
  ! This routine handles the potential failure of an individual Backward Euler step. Repeated calls
  ! to the BE step integrator are made with successivley smaller timesteps until success is achieved
  ! or ktsmx trials fail. The abundances are evolved using a Newton-Raphson iteration scheme to
  ! solve the equation (yt(i)-y(i))/tdel = ydot(i), where y(i) is the abundance at the beginning of
  ! the iteration, yt(i) is the trial abundance at the end of the timestep, ydot(i) is the time
  ! derivative of the trial abundance calculated from reaction rates and tdel is the timestep.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: y, ydot, yo, yt
  Use conditions
  Use controls, Only: idiag, iheat, iweak, iweak0, kitmx, kmon, ktot, lun_diag, lun_stdout, &
    & tdel_maxmult, nzbatchmx, szbatch
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_tstep
  Use xnet_interface, Only: step_be, t9rhofind, update_iweak
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep

  ! Input/Output variables
  Integer, Intent(inout) :: its(:) ! On input,  = 0 indicates active zone
                                   ! On output, = 1 if zone fails to converge

  ! Local variables
  Integer, Parameter :: ktsmx = 10
  Integer :: kts, j, izb, izone
  Integer :: inr(nzbatchmx)
  Integer :: mykts(nzbatchmx)

  start_timer = xnet_wtime()
  timer_tstep = timer_tstep - start_timer

  ! If the zone has previously converged or failed, do not iterate
  Where ( its /= 0 )
    inr = -1
  ElseWhere
    inr = 0
  EndWhere

  !-------------------------------------------------------------------------------------------------
  ! For each trial timestep, tdel, the NR iteration is attempted.
  ! The possible results for a timestep are a converged set of yt or a failure to converge.
  ! If convergence fails, retry with the trial timestep reduced by tdel_maxmult.
  !-------------------------------------------------------------------------------------------------
  mykts = 0
  Do kts = 1, ktsmx

    ! Attempt Backward Euler integration over desired timestep
    Call step_be(kstep,inr)

    ! Record number of NR iterations
    Where ( its == 0 )
      Where ( inr > 0 )
        kmon(2,:) = inr
        ktot(2,:) = ktot(2,:) + inr
      ElseWhere ( inr == 0 )
        kmon(2,:) = kitmx
        ktot(2,:) = ktot(2,:) + kitmx
      EndWhere
    EndWhere

    Do izb = 1, nzbatchmx
      ! If integration fails, reset abundances, reduce timestep and retry.
      If ( inr(izb) == 0 ) Then
        tdel(izb) = tdel(izb) / tdel_maxmult
        tt(izb) = t(izb) + tdel(izb)
        yet(izb) = ye(izb)
        yt(:,izb) = y(:,izb)
        iweak(izb) = iweak0

      ! If integration is successfull, flag the zone for removal from zone loop
      ElseIf ( inr(izb) > 0 ) Then
        mykts(izb) = kts
      EndIf
    EndDo

    ! Reset temperature and density for failed integrations
    Call t9rhofind(kstep,tt,ntt,t9t,rhot,mask_in = (inr == 0))
    If ( iheat > 0 ) Then
      Where ( inr == 0 )
        t9t = t9
      EndWhere
    EndIf
    Call update_iweak(t9t,mask_in = (inr == 0))

    ! Log the failed integration attempts
    If ( idiag >= 1 ) Then
      Do izb = 1, nzbatchmx
        If ( inr(izb) == 0 ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a,i5,i3,3es12.4)") 'BE TS Reduce',izone,kts,tt(izb),tdel(izb),tdel_maxmult
        EndIf
      EndDo
    EndIf

    If ( all( inr /= 0 ) ) Then
      Exit
    EndIf
  EndDo

  ! Mark TS convergence only for zones which haven't previously failed or converged
  Do izb = 1, nzbatchmx
    If ( inr(izb) > 0 ) Then
      nto(izb) = nt(izb)   ; nt(izb) = ntt(izb)
      to(izb) = t(izb)     ; t(izb) = tt(izb)
      t9o(izb) = t9(izb)   ; t9(izb) = t9t(izb)
      rhoo(izb) = rho(izb) ; rho(izb) = rhot(izb)
      yeo(izb) = ye(izb)   ; ye(izb) = yet(izb)
      yo(:,izb) = y(:,izb) ; y(:,izb) = yt(:,izb)
      tdel_next(izb) = tdel(izb) * tdel_maxmult
    ElseIf ( inr(izb) == 0 ) Then
      its(izb) = 1
    EndIf
  EndDo

  ! Log TS success/failure
  Do izb = 1, nzbatchmx
    izone = izb + szbatch - 1
    If ( inr(izb) > 0 ) Then
      If ( idiag >= 1 ) Write(lun_diag,"(a,2i5,2i3)") 'TS Success',kstep,izone,mykts(izb),inr(izb)
    ElseIf ( inr(izb) == 0 ) Then
      Write(lun_stdout,"(2i5,4es12.4,2i3)") kstep,izone,t(izb),tdel(izb),t9t(izb),rhot(izb),inr(izb),mykts(izb)
      Write(lun_stdout,*) 'Timestep retrys fail after ',mykts(izb),' attempts'
    EndIf
  EndDo

  kmon(1,:) = mykts
  ktot(1,:) = ktot(1,:) + mykts

  stop_timer = xnet_wtime()
  timer_tstep = timer_tstep + stop_timer

  Return
End Subroutine solve_be

Subroutine step_be(kstep,inr)
  !-------------------------------------------------------------------------------------------------
  ! This routine attempts to integrate a single Backward Euler step for the timestep tdel.
  ! If successful, inr = 1
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: y, ydot, yt
  Use conditions, Only: cv, rhot, t9, t9dot, t9t, tdel
  Use controls, Only: iconvc, idiag, iheat, ijac, kitmx, lun_diag, tolc, tolm, tolt9, ymin, &
    & nzbatchmx, szbatch
  Use nuc_number, Only: ny
  Use nuclear_data, Only: aa, nname
  Use thermo_data, Only: nh
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_nraph
  Use xnet_eos, Only: eos_cv
  Use xnet_interface, Only: cross_sect, jacobian_bksub, jacobian_build, jacobian_decomp, &
    & jacobian_solve, yderiv, yderiv_and_jacobian_build
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep

  ! Input/Output variables
  Integer, Intent(inout) :: inr(:) ! On input,  = 0 indicates active zone;
                                   !            =-1 indicates inactive zone
                                   ! On output, > 0 indicates # NR iterations if converged

  ! Local variables
  Integer :: irdymx, idymx
  Integer :: i, j, k, kit, izb, izone
  Real(dp) :: testc, testc2, testm, testn(nzbatchmx), toln(nzbatchmx)
  Real(dp) :: xtot(nzbatchmx), xtot_init(nzbatchmx), rdt(nzbatchmx)
  Real(dp) :: yrhs(ny,nzbatchmx), dy(ny,nzbatchmx), reldy(ny)
  Real(dp) :: t9rhs(nzbatchmx), dt9(nzbatchmx), relt9, mult
  Logical :: lznr(nzbatchmx)

  start_timer = xnet_wtime()
  timer_nraph = timer_nraph - start_timer

  ! Create mask for active zone integrations
  lznr = ( inr == 0 )

  ! Calculate the thermodynamic factors necessary for reaction rates, including
  ! screening and the reaction rates, if thermodynamic conditions are changing.
  If ( iheat == 0 ) Call cross_sect(mask_in = (nh > 1 .and. lznr))

  ! Calculate initial total mass fraction
  Do izb = 1, nzbatchmx
    If ( lznr(izb) ) Then
      xtot_init(izb) = sum(aa*y(:,izb)) - 1.0
      rdt(izb) = 1.0 / tdel(izb)
    EndIf
  EndDo

  ! The Newton-Raphson iteration occurs for at most kitmx iterations.
  Do kit = 1, kitmx

    ! If thermodynamic conditions are changing, rates need to be udpated
    If ( iheat > 0 ) Then
      Call cross_sect(mask_in = lznr)
      Do izb = 1, nzbatchmx
        If ( lznr(izb) ) Then
          Call eos_cv(rhot(izb),t9t(izb),yt(:,izb),cv(izb))
        EndIf
      EndDo
    EndIf

    ! Calculate the reaction rates and abundance time derivatives
    ! Rebuild and update LU factorization of the jacobian every ijac iterations
    If ( mod(kit-1,ijac) == 0 ) Then
      Call yderiv_and_jacobian_build(rdt,-1.0d0,mask_in = lznr)
      Call jacobian_decomp(kstep,mask_in = lznr)
    Else
      Call yderiv(mask_in = lznr)
    EndIf

    ! Calculate equation to zero
    Do izb = 1, nzbatchmx
      If ( lznr(izb) ) Then
        yrhs(:,izb) = (y(:,izb)-yt(:,izb))*rdt(izb) + ydot(:,izb)
        If ( iheat > 0 ) t9rhs(izb) = (t9(izb)-t9t(izb))*rdt(izb) + t9dot(izb)
      EndIf
    EndDo
    If ( idiag >= 2 ) Then
      Do izb = 1, nzbatchmx
        If ( lznr(izb) ) Then
          izone = izb + szbatch - 1
          Write(lun_diag,"(a3,2i5,es14.7)") 'RHS',kstep,izone,rdt(izb)
          Write(lun_diag,"(a5,4es23.15)") (nname(i),yrhs(i,izb),ydot(i,izb),yt(i,izb),y(i,izb),i=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(a5,4es23.15)") '   T9',t9rhs(izb),t9dot(izb),t9t(izb),t9(izb)
        EndIf
      EndDo
    EndIf

    ! Solve the jacobian and calculate the changes in abundances, dy
!   Call jacobian_solve(kstep,yrhs,dy,t9rhs,dt9,mask_in = lznr)
!   Call jacobian_decomp(kstep,mask_in = lznr)
    Call jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in = lznr)

    ! Evolve the abundances and calculate convergence tests
    Do izb = 1, nzbatchmx
      If ( lznr(izb) ) Then
        yt(:,izb) = yt(:,izb) + dy(:,izb)
        Where ( yt(:,izb) < ymin )
          yt(:,izb) = 0.0
          reldy = 0.0
        ElseWhere
          reldy = abs(dy(:,izb) / yt(:,izb))
        EndWhere
        If ( iheat > 0 ) Then
          t9t(izb) = t9t(izb) + dt9(izb)
          If ( t9t(izb) /= 0.0 ) Then
            relt9 = abs(dt9(izb) / t9t(izb))
          Else
            relt9 = 0.0
          EndIf
        Else
          relt9 = 0.0
        EndIf
        If ( idiag >= 3 ) Then
          izone = izb + szbatch - 1
          irdymx = maxloc(reldy,dim=1)
          idymx = maxloc(dy(:,izb),dim=1)
          Write(lun_diag,"(a3,2i5,i3,2(a5,2es23.15))") &
            & ' dY',kstep,izone,kit,nname(idymx),dy(idymx,izb),y(idymx,izb),nname(irdymx),reldy(irdymx),y(irdymx,izb)
          If ( idiag >= 4 ) Write(lun_diag,"(a5,5es12.4)") &
            & (nname(k),yt(k,izb),dy(k,izb),reldy(k),(aa(k)*dy(k,izb)),(aa(k)*yt(k,izb)),k=1,ny)
          If ( iheat > 0 ) Write(lun_diag,"(a3,2i5,i3,5x,2es23.15,5x,es12.4)") &
            & 'dT9',kstep,izone,kit,dt9(izb),t9t(izb),relt9
        EndIf

        !-------------------------------------------------------------------------------------------
        ! There are 3 included convergence tests:
        ! testc which measures relative changes,
        ! testc2 which measures total abundance changes, and
        ! testm which tests mass conservation.
        !-------------------------------------------------------------------------------------------
        testc  = sum(reldy)
        testc2 = sum(aa*dy(:,izb))
        xtot   = sum(aa*yt(:,izb)) - 1.0
        testm  = xtot(izb) - xtot_init(izb)
        If ( idiag >= 2 ) Write(lun_diag,"(a,3i5,3es14.6)") 'NR',kstep,izone,kit,testm,testc,testc2

        !-------------------------------------------------------------------------------------------
        ! testc is the most stringent test, and requires the most iterations.
        ! testm is the most lax, and therefore the fastest, often requiring only one iteration.
        ! Considering the uncertainties of the reaction rates, it is doubtful that the
        ! increased precision of testc is truly increased accuracy.
        !-------------------------------------------------------------------------------------------
        ! Ordinarily, test for true convergence
        If ( iconvc /= 0 ) Then
          testn(izb) = testc
          toln(izb) = tolc

        ! Otherwise, use mass conservation for convergence condition
        Else
          testn(izb) = testm
          toln(izb) = tolm
        EndIf

        ! If converged, exit NR loop
        If ( abs(testn(izb)) <= toln(izb) .and. relt9 < tolt9 ) Then
          inr(izb) = kit
          lznr(izb) = .false.
        EndIf
      EndIf
    EndDo

    ! Check that all zones are converged
    If ( all( inr /= 0 ) ) Exit
  EndDo

  If ( idiag >= 1 ) Then
    Do izb = 1, nzbatchmx
      izone = izb + szbatch - 1
      If ( inr(izb) > 0 ) Then
        Write(lun_diag,"(a,3i5,3es12.4)") 'Conv',kstep,izone,inr(izb),xtot(izb),testn(izb),toln(izb)
      ElseIf ( inr(izb) == 0 ) Then
        Write(lun_diag,"(a,3i5,2es10.2)") 'BE Failure',izone,inr(izb),kitmx,xtot(izb),testn(izb)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_nraph = timer_nraph + stop_timer

  Return
End Subroutine step_be