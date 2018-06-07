!***************************************************************************************************
! common.f90 10/18/17
! This file contains subroutines shared by the various integrators.
!***************************************************************************************************

Subroutine timestep(kstep,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine calculates the trial timestep.
  ! For tdel >0, this calculation is based on the relative changes of the abundances during the
  !              previous step in the evolution and the integrators estimate of the next timestep.
  !          =0, the timestep is calculated using the time derivatives of the abundances based on
  !              reaction rates.
  !          <0, the timestep is held constant, tdel=-tdel_old.
  ! There is also the provision for limiting the timestep in the event that the thermodynamic
  ! conditions are changing too rapidly.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: y, yo, yt, ydot
  Use conditions, Only: t, tt, tdel, tdel_next, tdel_old, t9, t9o, t9t, rho, rhot, &
    & t9dot, cv, nt, ntt, ints, intso
  Use controls, Only: changemx, changemxt, idiag, iheat, iweak, iweak0, lun_diag, yacc, &
    & nzbatchmx, szbatch, lzactive, iscrn
  Use nuc_number, Only: ny
  Use nuclear_data, Only: nname
  Use thermo_data, Only: tstop, tdelstart, nh, th, t9h, rhoh
  Use xnet_eos, Only: eos_cv
  Use xnet_interface, Only: cross_sect, t9rhofind, update_iweak, yderiv
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Integer, Intent(in) :: kstep

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i, j, n, nnew, nold, izb, izone
  Real(dp) :: ydotoy(ny), t9dotot9
  Real(dp) :: changest, changeth, dt, dtherm(nzbatchmx), ttest(nzbatchmx)
  Real(dp) :: tdel_stop, tdel_deriv, tdel_fine, tdel_t9, tdel_init
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  ! Retain old values of timestep and thermo and calculate remaining time
  changeth = 0.1
  changest = 0.01
  Where ( mask )
    tdel_old = tdel
  EndWhere
  Call cross_sect(mask_in = ( mask .and. tdel_old == 0.0 ))
  If ( iheat > 0 .and. iscrn <=0 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) .and. tdel_old(izb) == 0.0 ) Then
        Call eos_cv(rhot(izb),t9t(izb),yt(:,izb),cv(izb))
      EndIf
    EndDo
  EndIf
  Call yderiv(mask_in = ( mask .and. tdel_old == 0.0 ))

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      izone = izb + szbatch - 1

      tdel_stop = tstop(izb) - t(izb)
      tdel_fine = 0.0
      tdel_t9 = 0.0
      tdel_init = min(changest*tdelstart(izb),tdel_stop)

      ! If this is not the initial timestep, calculate timestep from changes in last timestep.
      If ( tdel_old(izb) > 0.0 ) Then
        Where ( y(:,izb) > yacc )
          ydotoy = abs((y(:,izb)-yo(:,izb))/y(:,izb))
        ElseWhere
          ydotoy = 0.0
        EndWhere
        ints(izb) = maxloc(ydotoy,dim=1)
        If ( ydotoy(ints(izb)) /= 0.0 ) Then
          tdel_deriv = tdel_old(izb) * changemx/ydotoy(ints(izb))
        Else
          tdel_deriv = tdel_next(izb)
        EndIf
        If ( iheat > 0 ) Then
          t9dotot9 = abs((t9(izb)-t9o(izb))/t9(izb))
          If ( t9dotot9 > 0.0 ) Then
            tdel_t9 = tdel_old(izb) * changemxt/t9dotot9
          Else
            tdel_t9 = tdel_next(izb)
          EndIf
        Else
          tdel_t9 = tdel_next(izb)
        EndIf
        tdel(izb) = min( tdel_deriv, tdel_stop, tdel_next(izb), tdel_t9 )

      ! If this is an initial timestep, yo does not exist, so calculate timestep directly from derivatives.
      ElseIf ( tdel_old(izb) == 0.0 ) Then
        If ( nh(izb) > 2 ) tdel_stop = 1.0e-4*tdel_stop
        intso(izb) = 0

        Where ( y(:,izb) > yacc )
          ydotoy = abs(ydot(:,izb)/y(:,izb))
        ElseWhere
          ydotoy = 0.0
        EndWhere

        ! If derivatives are non-zero, Use Y/(dY/dt).
        ints(izb) = maxloc(ydotoy,dim=1)
        If ( ydotoy(ints(izb)) /= 0.0 ) Then
          tdel_deriv = changemx/ydotoy(ints(izb))

          ! For unevolved initial abundances, as often found in test problems,
          ! tdel_deriv may produce large jumps from zero abundance in the first
          ! timestep. While generally inconsequential, tdel_fine limits these
          ! to the accuracy abundance limit.
          tdel_fine = changest*tdel_deriv

          ! If derivatives are zero, take a small step.
        Else
          tdel_deriv = tdel_stop
          tdel_fine = tdel_stop
        EndIf
        If ( iheat > 0 ) Then
          tdel_t9 = changemxt * abs(t9(izb)/t9dot(izb))
        Else
          tdel_t9 = tdel_stop
        EndIf
        tdel(izb) = min( tdel_stop, tdel_deriv, tdel_fine, tdel_t9 )

        ! Use the user-defined initial timestep if it is larger than the estimated timestep
        tdel(izb) = max( tdel_init, tdel(izb) )

        ! Keep timestep constant
      Else
        tdel(izb) = -tdel_old(izb)
      EndIf

      If ( idiag >= 1 ) Write(lun_diag,"(a4,2i5,7es12.4,i5)") &
        & 'tdel',kstep,izone,tdel(izb),tdel_deriv,tdel_old(izb),tdel_stop,tdel_next(izb),tdel_fine,tdel_t9,ints(izb)
!     If ( idiag >= 2 ) Write(lun_diag,"(a5,i4,2es12.4)") (nname(k),k,y(k,izb),ydotoy(k),k=1,ny)

      ! Retain the index of the species setting the timestep
      If ( ints(izb) /= intso(izb) ) Then
        If ( idiag >= 1 ) Write(lun_diag,"(a4,a5,3es23.15)") 'ITC ',nname(ints(izb)),y(ints(izb),izb),t(izb),tdel(izb)
        intso(izb) = ints(izb)
      EndIf
    EndIf
  EndDo

  ! For varying temperature and density, capture thermodynamic features
  If ( iheat == 0 ) Then

    ! Make sure to not skip any features in the temperature or density profiles by checking
    ! for profile monotonicity between t and t+del
    Where ( mask .and. nh > 1 )
      ttest = min(t + tdel, tstop)
    EndWhere
    Call t9rhofind(kstep,ttest,ntt,t9t,rhot,mask_in = ( mask .and. nh > 1 ))
    Do izb = 1, nzbatchmx
      If ( mask(izb) .and. nh(izb) > 1 .and. ntt(izb)-1 > nt(izb) ) Then
        Do j = nt(izb), ntt(izb)-1
          If ( t9h(j,izb) > t9(izb) .and. t9h(j,izb) > t9t(izb) ) Then
            tdel(izb) = th(j,izb) - t(izb)
            Exit
          ElseIf ( t9h(j,izb) < t9(izb) .and. t9h(j,izb) < t9t(izb) ) Then
            tdel(izb) = th(j,izb) - t(izb)
            Exit
          ElseIf ( rhoh(j,izb) > rho(izb) .and. rhoh(j,izb) > rhot(izb) ) Then
            tdel(izb) = th(j,izb) - t(izb)
            Exit
          ElseIf ( rhoh(j,izb) < rho(izb) .and. rhoh(j,izb) < rhot(izb) ) Then
            tdel(izb) = th(j,izb) - t(izb)
            Exit
          EndIf
        EndDo
      EndIf
    EndDo

    ! Limit timestep if fractional density change is larger than changeth (10% by default)
    ! or fraction temperature change is larger than 0.1*changeth (1% by default)
    dtherm = 0.0
    Do i = 1, 10
      Where ( mask .and. nh > 1 )
        ttest = min(t + tdel, tstop)
      EndWhere
      Call t9rhofind(kstep,ttest,ntt,t9t,rhot,mask_in = ( mask .and. nh > 1 ))
      Where ( mask .and. nh > 1 .and. t9 > 0.0 )
        dtherm = 10.0*abs(t9t-t9)/t9 + abs(rhot-rho)/rho
      EndWhere
      If ( all( dtherm < changeth ) ) Exit
      tdel(izb) = 0.5*tdel(izb)
    EndDo
    Do izb = 1, nzbatchmx
      If ( mask(izb) .and. nh(izb) > 1 ) Then
        izone = izb + szbatch - 1
        If ( dtherm(izb) >= changeth ) Write(lun_diag,"(a,i2,a,i5,3es23.15)") &
          & 'Error in Thermo variations after ',i,'reductions',izone,tdel(izb),t9t(izb),rhot(izb)
        If ( idiag >= 1 ) Write(lun_diag,"(a5,2i5,2es23.15)") 'T9del',kstep,izone,tdel(izb),dtherm(izb)
      EndIf
    EndDo
  EndIf

  Where ( mask )
    tt = min(t + tdel, tstop)
    iweak = iweak0
  EndWhere

  ! Check for any changes to iweak
  Call update_iweak(t9t,mask_in = mask)

  Return
End Subroutine timestep

Subroutine partf(t9,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine calculates the nuclear partition functions as a function of temperature.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: idiag, iheat, lun_diag, nzbatchmx, lzactive
  Use nuc_number, Only: ny
  Use part_funct_data, Only: dlngdt9, g, gg, ng, t9i
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Real(dp), Intent(in) :: t9(:)

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i, ii, izb
  Real(dp) :: rdt9
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      Do i = 1, ng
        If ( t9(izb) <= t9i(i) ) Exit
      EndDo
      ii = i

      ! Linear interpolation in log-space
      Select Case (ii)
      Case (1)
        gg(1:ny,izb) = g(1,1:ny)
      Case (ng+1)
        gg(1:ny,izb) = g(ng,1:ny)
      Case Default
        rdt9 = (t9(izb)-t9i(ii-1)) / (t9i(ii)-t9i(ii-1))
        gg(1:ny,izb) = exp( rdt9*log(g(ii,1:ny)) + (1.0-rdt9)*log(g(ii-1,1:ny)) )
      End Select
      gg(0,izb) = 1.0 ! placeholder for non-nuclei, gamma-rays, etc.

      If ( iheat > 0 ) Then
        Select Case (ii)
        Case (1)
          dlngdt9(1:ny,izb) = log(g(2,1:ny)/g(1,1:ny)) / (t9i(2)-t9i(1))
        Case (ng+1)
          dlngdt9(1:ny,izb) = log(g(ng,1:ny)/g(ng-1,1:ny)) / (t9i(ng)-t9i(ng-1))
        Case Default
          dlngdt9(1:ny,izb) = log(g(ii,1:ny)/g(ii-1,1:ny)) / (t9i(ii)-t9i(ii-1))
        End Select
        dlngdt9(0,izb) = 0.0
      EndIf
    EndIf
  EndDo

! If ( idiag >= 1 ) Then
!   Do izb = 1, nzbatchmx
!     If ( mask(izb) ) Then
!       Write(lun_diag,"(a5,i3,es14.7)") 'PartF',ii,t9(izb)
!       Write(lun_diag,"(5(i4,es12.4))") (i, gg(i,izb), i=1,ny)
!     EndIf
!   EndDo
! EndIf

  Return
End Subroutine partf

Subroutine update_iweak(t9,mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine updates the value of the iweak flag to control the treatment of strong reactions.
  !-------------------------------------------------------------------------------------------------
  Use controls, Only: iweak0, iweak, lun_stdout, t9min, nzbatchmx, lzactive
  Use xnet_types, Only: dp
  Implicit None

  ! Input variables
  Real(dp), Intent(in) :: t9(:)

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: izb
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then

      ! Turn off strong interactions if the temperature is less than t9min
      If ( t9(izb) < t9min ) Then
        If ( iweak(izb) /= 0 ) Then
          iweak(izb) = -1
        Else
          iweak(izb) = iweak0
          Write(lun_stdout,*) 'Warning: Strong reactions ignored for T9 < T9min.'
        EndIf
      EndIf
    EndIf
  EndDo

  Return
End Subroutine update_iweak

Subroutine yderiv(mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine calculates time derivatives for each nuclear species. This calculation is performed
  ! by looping over nuclei and summing the reaction rates for each reaction which effects that
  ! nucleus (see Equation 10 in Hix & Meyer (2006)).
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: yt, ydot
  Use conditions, Only: cv, t9t, t9dot
  Use controls, Only: idiag, iheat, lun_diag, nzbatchmx, szbatch, lzactive
  Use cross_sect_data, Only: csect1, csect2, csect3, n1i, n2i, n3i
  Use nuc_number, Only: ny
  Use nuclear_data, Only: nname, mex
  Use reac_rate_data, Only: a1, a2, a3, b1, b2, b3, la, le, mu1, mu2, mu3, &
    & n11, n21, n22, n31, n32, n33
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_deriv
  Use xnet_types, Only: dp
  Implicit None

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer :: i, j, i0, i1, la1, le1, la2, le2, la3, le3, izb, izone
  Real(dp) :: s1, s2, s3, s11, s22, s33
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf

  start_timer = xnet_wtime()
  timer_deriv = timer_deriv - start_timer

  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then

      ! From the cross sections and the counting array, calculate the reaction rates
      b1(:,izb) = a1*csect1(mu1,izb)
      b2(:,izb) = a2*csect2(mu2,izb)
      b3(:,izb) = a3*csect3(mu3,izb)

      ! Calculate Ydot for each nucleus, summing over the reactions which affect it.
      Do i0 = 1, ny
        la1 = la(1,i0) ; la2 = la(2,i0) ; la3 = la(3,i0)
        le1 = le(1,i0) ; le2 = le(2,i0) ; le3 = le(3,i0)

        ! Sum over the reactions with 1 reactant
        s1 = 0.0
        Do i1 = la1, le1
          s11 = b1(i1,izb)*yt(n11(i1),izb)
          s1 = s1 + s11
        EndDo

        ! Sum over the reactions with 2 reactants
        s2 = 0.0
        Do i1 = la2, le2
          s22 = b2(i1,izb)*yt(n21(i1),izb)*yt(n22(i1),izb)
          s2 = s2 + s22
        EndDo

        ! Sum over the reactions with 3 reactants
        s3 = 0.0
        Do i1 = la3, le3
          s33 = b3(i1,izb)*yt(n31(i1),izb)*yt(n32(i1),izb)*yt(n33(i1),izb)
          s3 = s3 + s33
        EndDo

        ! Sum the 3 components of Ydot
        ydot(i0,izb) = s1 + s2 + s3
      EndDo

      If ( iheat > 0 ) Then
        ! Surprisingly, this seems to perform better than the DGEMV below
        s1 = 0.0
        Do i0 = 1, ny
          s1 = s1 + mex(i0)*ydot(i0,izb)
        EndDo
        t9dot(izb) = -s1 / cv(izb)
      EndIf
    EndIf
  EndDo
! If ( iheat > 0 ) Then
!   Call dgemv('T',ny,nzbatchmx,1.0,ydot,ny,mex,1,0.0,t9dot,1)
!   Where ( mask )
!     t9dot = -t9dot / cv
!   EndWhere
! EndIf

  ! Separate loop for diagnostics so compiler can properly vectorize
  If ( idiag >= 3 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a,i5)") 'YDERIV',izone
        Do i0 = 1, ny
          la1 = la(1,i0) ; la2 = la(2,i0) ; la3 = la(3,i0)
          le1 = le(1,i0) ; le2 = le(2,i0) ; le3 = le(3,i0)
          s1 = 0.0 ; s2 = 0.0 ; s3 = 0.0
          Do i1 = la1, le1
            s11 = b1(i1,izb)*yt(n11(i1),izb)
            s1 = s1 + s11
          EndDo
          Do i1 = la2, le2
            s22 = b2(i1,izb)*yt(n21(i1),izb)*yt(n22(i1),izb)
            s2 = s2 + s22
          EndDo
          Do i1 = la3, le3
            s33 = b3(i1,izb)*yt(n31(i1),izb)*yt(n32(i1),izb)*yt(n33(i1),izb)
            s3 = s3 + s33
          EndDo
          If ( idiag >= 5 ) Then
            Write(lun_diag,"(a3,a6,i4)") 'NUC',nname(i0),i0
            Write(lun_diag,"(3a5,'  1  ',4es23.15)") &
              & ((nname(n1i(j,mu1(i1))),j=1,3),b1(i1,izb),yt(n11(i1),izb), &
              & b1(i1,izb)*yt(n11(i1),izb),a1(i1),i1=la1,le1)
            Write(lun_diag,*) '1->',nname(i0),la1,le1,s1
            Write(lun_diag,"(4a5,4es23.15)") &
              & ((nname(n2i(i,mu2(i1))),i=1,4),b2(i1,izb),yt(n21(i1),izb),yt(n22(i1),izb), &
              & b2(i1,izb)*yt(n21(i1),izb)*yt(n22(i1),izb),i1=la2,le2)
            Write(lun_diag,*) '2->',nname(i0),la2,le2,s2
            Write(lun_diag,"(3a5,'  3  ',5es23.15)") &
              & ((nname(n3i(i,mu3(i1))),i=1,3),b3(i1,izb),yt(n31(i1),izb),yt(n32(i1),izb),yt(n33(i1),izb), &
              & b3(i1,izb)*yt(n31(i1),izb)*yt(n32(i1),izb)*yt(n33(i1),izb),i1=la3,le3)
            Write(lun_diag,*) '3->',nname(i0),la3,le3,s3
          EndIf
          Write(lun_diag,"(a4,a5,2es23.15,3es11.3)") 'YDOT',nname(i0),yt(i0,izb),ydot(i0,izb),s1,s2,s3
        EndDo
        If ( iheat > 0 ) Write(lun_diag,"(5x,a5,2es23.15)") '   T9',t9t(izb),t9dot(izb)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_deriv = timer_deriv + stop_timer

  Return
End Subroutine yderiv

Subroutine cross_sect(mask_in)
  !-------------------------------------------------------------------------------------------------
  ! This routine calculates the cross section for each reaction.
  !-------------------------------------------------------------------------------------------------
  Use abundances, Only: yt
  Use conditions, Only: t, rhot, t9t, yet
  Use controls, Only: idiag, iheat, iscrn, iweak, lun_diag, nzbatchmx, szbatch, lzactive
  Use cross_sect_data
  Use nuc_number, Only: ny
  Use nuclear_data, Only: nname
  Use part_funct_data, Only: gg, dlngdt9
  Use screening_data, Only: h1, h2, h3, dh1dt9, dh2dt9, dh3dt9
  Use timers, Only: xnet_wtime, start_timer, stop_timer, timer_scrn, timer_csect
  Use xnet_eos, Only: y_moment
  Use xnet_interface, Only: ffn_rate, partf, screening
  Use xnet_types, Only: dp
  Implicit None

  ! Optional variables
  Logical, Optional, Target, Intent(in) :: mask_in(:)

  ! Local variables
  Integer, Parameter :: dgemm_nzbatch = 2 ! Min number of zones to use dgemm, otherwise use dgemv
  Real(dp) :: t09(7,nzbatchmx), dt09(7,nzbatchmx)
  Real(dp) :: ene(nzbatchmx), ytot, abar, zbar, z2bar, zibar
  Real(dp) :: rffn(nffn,nzbatchmx)        ! FFN reaction rates
  Real(dp) :: dlnrffndt9(nffn,nzbatchmx)  ! log FFN reaction rates derivatives
! Real(dp) :: rnnu(nnnu,4,nzbatchmx)      ! Neutrino reaction rates                             !NNU
  Real(dp) :: rhot2, rhot3
  Integer :: j, k, izb, izone, nzmask
  Logical, Pointer :: mask(:)

  If ( present(mask_in) ) Then
    mask => mask_in
  Else
    mask => lzactive
  EndIf
  nzmask = count(mask)

  start_timer = xnet_wtime()
  timer_scrn = timer_scrn - start_timer

  ! Calculate the screening terms (also updates electron fraction to be consistent with yt)
  If ( iscrn >= 1 ) Then
    Call screening(mask_in = mask)
  Else
    h1 = 0.0
    h2 = 0.0
    h3 = 0.0
!   h4 = 0.0
    If ( iheat > 0 ) Then
      dh1dt9 = 0.0
      dh2dt9 = 0.0
      dh3dt9 = 0.0
!     dh4dt9 = 0.0
    EndIf
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        Call y_moment(yt(:,izb),yet(izb),ytot,abar,zbar,z2bar,zibar)
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_scrn = timer_scrn + stop_timer

  start_timer = xnet_wtime()
  timer_csect = timer_csect - start_timer

  ! Calculate necessary thermodynamic moments
  t09 = 0.0
  Where ( mask )
    ene = yet*rhot
    t09(1,:) = 1.0
    t09(2,:) = t9t**(-1)
    t09(3,:) = t9t**(-1.0/3.0)
    t09(4,:) = t9t**(+1.0/3.0)
    t09(5,:) = t9t
    t09(6,:) = t9t**(+5.0/3.0)
    t09(7,:) = log(t9t)
  EndWhere

  ! Calculate the REACLIB exponent polynomials, adding screening terms
  If ( nzmask >= dgemm_nzbatch ) Then
    Call dgemm('T','N',nreac(1),nzbatchmx,7,1.0,rc1,7,t09,7,1.0,h1,nreac(1))
    Call dgemm('T','N',nreac(2),nzbatchmx,7,1.0,rc2,7,t09,7,1.0,h2,nreac(2))
    Call dgemm('T','N',nreac(3),nzbatchmx,7,1.0,rc3,7,t09,7,1.0,h3,nreac(3))
  Else
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        Call dgemv('T',7,nreac(1),1.0,rc1,7,t09(:,izb),1,1.0,h1(:,izb),1)
        Call dgemv('T',7,nreac(2),1.0,rc2,7,t09(:,izb),1,1.0,h2(:,izb),1)
        Call dgemv('T',7,nreac(3),1.0,rc3,7,t09(:,izb),1,1.0,h3(:,izb),1)
      EndIf
    EndDo
  EndIf

  If ( idiag >= 3 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a3,i5,4es23.15)") 'THR',izone,t9t(izb),rhot(izb),ene(izb),yet(izb)
      EndIf
    EndDo
  EndIf

  ! Calculate partition functions for each nucleus at t9t
  Call partf(t9t,mask_in = mask)

  ! If there are any FFN reactions, calculate their rates
  If ( nffn > 0 ) Call ffn_rate(nffn,t9t,ene,rffn,dlnrffndt9,mask_in = mask)                    !FFN

  ! If there are any neutrino-nucleus reactions, calculate their rates
! If ( nnnu > 0 ) Call nnu_rate(t,rnnu,mask_in = mask)                                          !NNU

  ! Calculate the csect for reactions with 1 reactant
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      Where ( irev1 == 1 )
        rpf1 = gg(n1i(2,:),izb) * gg(n1i(3,:),izb) / gg(n1i(1,:),izb)
      ElseWhere
        rpf1 = 1.0
      EndWhere
      If ( iweak(izb) > 0 ) Then
        Where ( iwk1 == 0 .or. iwk1 == 4 )
          csect1(:,izb) = rpf1 * exp(h1(:,izb))
        ElseWhere ( iwk1 == 1 )
          csect1(:,izb) = rpf1 * exp(h1(:,izb)) * ene(izb)
        ElseWhere ( iwk1 == 2 .or. iwk1 == 3 )  ! FFN reaction
          csect1(:,izb) = rffn(iffn,izb)
!       ElseWhere ( iwk1 == 7 )             ! Electron neutrino capture                         !NNU
!         csect1(:,izb) = rpf1 * rnnu(innu,1,izb)                                               !NNU
!       ElseWhere ( iwk1 == 8 )             ! Electron anti-neutrino capture                    !NNU
!         csect1(:,izb) = rpf1 * rnnu(innu,2,izb)                                               !NNU
        EndWhere
      ElseIf ( iweak(izb) < 0 ) Then
        Where ( iwk1 == 0 )
          csect1(:,izb) = 0.0
        ElseWhere ( iwk1 == 1 )
          csect1(:,izb) = rpf1 * exp(h1(:,izb)) * ene(izb)
        ElseWhere ( iwk1 == 4 )
          csect1(:,izb) = rpf1 * exp(h1(:,izb))
        ElseWhere ( iwk1 == 2 .or. iwk1 == 3 )  ! FFN reaction
          csect1(:,izb) = rffn(iffn,izb)
!       ElseWhere ( iwk1 == 7 )             ! Electron neutrino capture                         !NNU
!         csect1(:,izb) = rpf1 * rnnu(innu,1,izb)                                               !NNU
!       ElseWhere ( iwk1 == 8 )             ! Electron anti-neutrino capture                    !NNU
!         csect1(:,izb) = rpf1 * rnnu(innu,2,izb)                                               !NNU
        EndWhere
      Else
        Where ( iwk1 == 0 )
          csect1(:,izb) = rpf1 * exp(h1(:,izb))
        ElseWhere
          csect1(:,izb) = 0.0
        EndWhere
      EndIf
    EndIf
  EndDo

  ! Calculate the csect for reactions with 2 reactants
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
      Where ( irev2 == 1 )
        rpf2 = gg(n2i(3,:),izb) * gg(n2i(4,:),izb) / ( gg(n2i(1,:),izb) * gg(n2i(2,:),izb) )
      ElseWhere
        rpf2 = 1.0
      EndWhere
      If ( iweak(izb) > 0 ) Then
        Where ( iwk2 == 1 )
          csect2(:,izb) = rhot(izb) * rpf2 * exp(h2(:,izb)) * ene(izb)
        ElseWhere
          csect2(:,izb) = rhot(izb) * rpf2 * exp(h2(:,izb))
        EndWhere
      ElseIf ( iweak(izb) < 0 ) Then
        Where ( iwk2 == 0 )
          csect2(:,izb) = 0.0
        ElseWhere ( iwk2 == 1 )
          csect2(:,izb) = rhot(izb) * rpf2 * exp(h2(:,izb)) * ene(izb)
        ElseWhere
          csect2(:,izb) = rhot(izb) * rpf2 * exp(h2(:,izb))
        EndWhere
      Else
        Where ( iwk2 == 0 )
          csect2(:,izb) = rhot(izb) * rpf2 * exp(h2(:,izb))
        ElseWhere
          csect2(:,izb) = 0.0
        EndWhere
      EndIf
    EndIf
  EndDo

  ! Calculate the csect for reactions with 3 reactants
  Do izb = 1, nzbatchmx
    If ( mask(izb) ) Then
!     Where ( irev3 == 1 )
!       rpf3 = gg(n3i(4,:),izb) * gg(n3i(5,:),izb) / ( gg(n3i(1,:),izb) * gg(n3i(2,:),izb) * gg(n3i(3,:),izb) )
!     ElseWhere
!       rpf3 = 1.0
!     EndWhere
      rhot2 = rhot(izb)**2
      If ( iweak(izb) > 0 ) Then
        Where ( iwk3 == 1 )
          csect3(:,izb) = rhot2 * exp(h3(:,izb)) * ene(izb)
        ElseWhere
          csect3(:,izb) = rhot2 * exp(h3(:,izb))
        EndWhere
      ElseIf ( iweak(izb) < 0 ) Then
        Where ( iwk3 == 0 )
          csect3(:,izb) = 0.0
        ElseWhere ( iwk3 == 1 )
          csect3(:,izb) = rhot2 * exp(h3(:,izb)) * ene(izb)
        ElseWhere
          csect3(:,izb) = rhot2 * exp(h3(:,izb))
        EndWhere
      Else
        Where ( iwk3 == 0 )
          csect3(:,izb) = rhot2 * exp(h3(:,izb))
        ElseWhere
          csect3(:,izb) = 0.0
        EndWhere
      EndIf
      Where ( csect3(:,izb) < 1.0e-20 )
        csect3(:,izb) = 0.0
      EndWhere
    EndIf
  EndDo

  ! Calculate the csect for reactions with 4 reactants
! Call dgemm('T','N',nreac(4),nzbatchmx,7,1.0,rc4,7,t09,7,1.0,h4,nreac(4))
! Do izb = 1, nzbatchmx
!   If ( mask(izb) ) Then
!     Where ( irev4 == 1 )
!       rpf4 = gg(n4i(5,:),izb) * gg(n4i(6,:),izb) / ( gg(n4i(1,:),izb) * gg(n4i(2,:),izb) * gg(n4i(3,:),izb) * gg(n4i(4,:),izb ) )
!     ElseWhere
!       rpf4 = 1.0
!     EndWhere
!     rhot3 = rhot(izb)**3
!     If ( iweak(izb) > 0 ) Then
!       Where ( iwk4 == 1 )
!         csect4(:,izb) = rhot3 * exp(h4(:,izb)) * ene(izb)
!       ElseWhere
!         csect4(:,izb) = rhot3 * exp(h4(:,izb))
!       EndWhere
!     ElseIf ( iweak(izb) < 0 ) Then
!       Where ( iwk4 == 0 )
!         csect4(:,izb) = 0.0
!       ElseWhere ( iwk4 == 1 )
!         csect4(:,izb) = rhot3 * exp(h4(:,izb)) * ene(izb)
!       ElseWhere
!         csect4(:,izb) = rhot3 * exp(h4(:,izb))
!       EndWhere
!     Else
!       Where ( iwk4 == 0 )
!         csect4(:,izb) = rhot3 * exp(h4(:,izb))
!       ElseWhere
!         csect4(:,izb) = 0.0
!       EndWhere
!     EndIf
!   EndIf
! EndDo

  If ( iheat > 0 ) Then

    ! Calculate reaction rate coefficient derivatives
    dt09 = 0.0
    Where ( mask )
      dt09(1,:) = 0.0
      dt09(2,:) = -t9t**(-2)
      dt09(3,:) = -t9t**(-4.0/3.0) / 3.0
      dt09(4,:) = +t9t**(-2.0/3.0) / 3.0
      dt09(5,:) = 1.0
      dt09(6,:) = +t9t**(+2.0/3.0) * 5.0/3.0
      dt09(7,:) = 1.0 / t9t
    EndWhere

    ! Calculate the derivatives of REACLIB exponent polynomials, adding screening terms
    If ( nzmask >= dgemm_nzbatch ) Then
      Call dgemm('T','N',nreac(1),nzbatchmx,7,1.0,rc1,7,dt09,7,1.0,dh1dt9,nreac(1))
      Call dgemm('T','N',nreac(2),nzbatchmx,7,1.0,rc2,7,dt09,7,1.0,dh2dt9,nreac(2))
      Call dgemm('T','N',nreac(3),nzbatchmx,7,1.0,rc3,7,dt09,7,1.0,dh3dt9,nreac(3))
    Else
      Do izb = 1, nzbatchmx
        If ( mask(izb) ) Then
          Call dgemv('T',7,nreac(1),1.0,rc1,7,dt09(:,izb),1,1.0,dh1dt9(:,izb),1)
          Call dgemv('T',7,nreac(2),1.0,rc2,7,dt09(:,izb),1,1.0,dh2dt9(:,izb),1)
          Call dgemv('T',7,nreac(3),1.0,rc3,7,dt09(:,izb),1,1.0,dh3dt9(:,izb),1)
        EndIf
      EndDo
    EndIf

    ! 1-reactant reactions
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        Where ( irev1 == 1 )
          dlnrpf1dt9 = dlngdt9(n1i(2,:),izb) + dlngdt9(n1i(3,:),izb) - dlngdt9(n1i(1,:),izb)
        ElseWhere
          dlnrpf1dt9 = 0.0
        EndWhere
        Where ( iwk1 == 2 .or. iwk1 == 3 )  ! FFN reaction
          dcsect1dt9(:,izb) = rffn(iffn,izb)*dlnrffndt9(iffn,izb)
!       ElseWhere ( iwk1 == 7 )             ! Electron neutrino capture                         !NNU
!         dcsect1dt9(:,izb) = rnnu(innu,1,izb)*dlnrpf1dt9                                       !NNU
!       ElseWhere ( iwk1 == 8 )             ! Electron anti-neutrino capture                    !NNU
!         dcsect1dt9(:,izb) = rnnu(innu,2,izb)*dlnrpf1dt9                                       !NNU
        ElseWhere
          dcsect1dt9(:,izb) = csect1(:,izb)*(dh1dt9(:,izb)+dlnrpf1dt9)
        EndWhere
      EndIf
    EndDo

    ! 2-reactant reactions
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        Where ( irev2 == 1 )
          dlnrpf2dt9 = dlngdt9(n2i(3,:),izb) + dlngdt9(n2i(4,:),izb) - dlngdt9(n2i(1,:),izb) - dlngdt9(n2i(2,:),izb)
        ElseWhere
          dlnrpf2dt9 = 0.0
        EndWhere
        dcsect2dt9(:,izb) = csect2(:,izb)*(dh2dt9(:,izb)+dlnrpf2dt9)
      EndIf
    EndDo

    ! 3-reactant reactions
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
!       Where ( irev3 == 1 )
!         dlnrpf3dt9 = + dlngdt9(n3i(4,:),izb) + dlngdt9(n3i(5,:),izb) &
!           &          - dlngdt9(n3i(1,:),izb) - dlngdt9(n3i(2,:),izb) - dlngdt9(n3i(3,:),izb) 
!       ElseWhere
!         dlnrpf3dt9 = 0.0
!       EndWhere
!       dcsect3dt9(:,izb) = csect3(:,izb)*(dh3dt9(:,izb)+dlnrpf3dt9)
        dcsect3dt9(:,izb) = csect3(:,izb)*dh3dt9(:,izb)
      EndIf
    EndDo

    ! 4-reactant reactions
!   Call dgemm('T','N',nreac(4),nzbatchmx,7,1.0,rc4,7,dt09,7,1.0,dh4dt9,nreac(4))
!   Do izb = 1, nzbatchmx
!     If ( mask(izb) ) Then
!       Where ( irev4 == 1 )
!         dlnrpf4dt9 = + dlngdt9(n4i(5,:),izb) + dlngdt9(n4i(6,:),izb) &
!           &          - dlngdt9(n4i(1,:),izb) - dlngdt9(n4i(2,:),izb) - dlngdt9(n4i(3,:),izb) - dlngdt9(n4i(4,:),izb)
!       ElseWhere
!         dlnrpf4dt9= 0.0
!       EndWhere
!       dcsect4dt9(:,izb) = csect4(:,izb)*(dh4dt9(:,izb)+dlnrpf4dt9)
!     EndIf
!   EndDo
  EndIf

  If ( idiag >= 5 ) Then
    Do izb = 1, nzbatchmx
      If ( mask(izb) ) Then
        izone = izb + szbatch - 1
        Write(lun_diag,"(a,i5)") 'CSect',izone
        Write(lun_diag,"(a,i5)") 'CSect1',nreac(1)
        Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
          & (k,nname(n1i(1,k)),'-->',(nname(n1i(j,k)),j=2,3),'+++', &
          & iwk1(k),irev1(k),ires1(k),h1(k,izb),csect1(k,izb), k=1,nreac(1))
        Write(lun_diag,"(a,i5)") 'CSect2',nreac(2)
        Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
          & (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,4), &
          & iwk2(k),irev2(k),ires2(k),h2(k,izb),csect2(k,izb),k=1,nreac(2))
        Write(lun_diag,"(a,i5)") 'CSect3',nreac(3)
        Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
          & (k,(nname(n3i(j,k)),j=1,3),'-->','+++', &
          & iwk3(k),irev3(k),ires3(k),h3(k,izb),csect3(k,izb),k=1,nreac(3))
!       Write(lun_diag,"(a,i5)") 'CSect4',nreac(4)
!       Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
!         & (k,(nname(n4i(j,k)),j=1,4),'-->','+++', &
!         & iwk4(k),irev4(k),ires4(k),h4(k,izb),csect4(k,izb),k=1,nreac(4))
        If ( iheat > 0 ) Then
          Write(lun_diag,"(a,i5)") 'dCSect1/dT9',nreac(1)
          Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
            & (k,nname(n1i(1,k)),'-->',(nname(n1i(j,k)),j=2,3),'+++',&
            & iwk1(k),irev1(k),ires1(k),dh1dt9(k,izb),dcsect1dt9(k,izb), k=1,nreac(1))
          Write(lun_diag,"(a,i5)") 'dCSect2/dT9',nreac(2)
          Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
            & (k,(nname(n2i(j,k)),j=1,2),'-->',(nname(n2i(j,k)),j=3,4), &
            & iwk2(k),irev2(k),ires2(k),dh2dt9(k,izb),dcsect2dt9(k,izb),k=1,nreac(2))
          Write(lun_diag,"(a,i5)") 'dCSect3/dT9',nreac(3)
          Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
            & (k,(nname(n3i(j,k)),j=1,3),'-->','+++', &
            & iwk3(k),irev3(k),ires3(k),dh3dt9(k,izb),dcsect3dt9(k,izb),k=1,nreac(3))
!         Write(lun_diag,"(a,i5)") 'dCSect4/dT9',nreac(4)
!         Write(lun_diag,"(i5,5a5,3i3,2es23.15)") &
!           & (k,(nname(n4i(j,k)),j=1,4),'-->','+++', &
!           & iwk4(k),irev4(k),ires4(k),dh4dt9(k,izb),dcsect4dt9(k,izb),k=1,nreac(4))
        EndIf
      EndIf
    EndDo
  EndIf

  stop_timer = xnet_wtime()
  timer_csect = timer_csect + stop_timer

  Return
End Subroutine cross_sect
