module actual_burner_module

  implicit none

contains

  subroutine actual_burner(state_in, state_out, dt, time)

    use eos_module
!    use Hydro_interface, ONLY : Hydro_detectShock
    use network
    use NSE_data
    use fundamental_constants_module
    use meth_params_module, only: small_temp

    implicit none

    type (eos_t) :: state_in, state_out
    double precision :: dt, time
    
    integer                    :: i, j, k, iref, l, n
    integer                    :: blockID, istat
    double precision           :: dx

    double precision :: shock
    double precision :: react_proximity

    double precision :: flamewidth
    integer :: ii, jj, kk, maxdi

    logical :: ignition_test, bcastMask, spark_block
    double precision :: spark_block_radius, radius
    logical, allocatable, dimension(:) :: ignition_conditions
    integer :: ignition_num
    double precision, allocatable, dimension(:) :: ignition_coords

    integer :: det_num
    double precision, allocatable, dimension(:) :: det_xCoord, det_yCoord, det_zCoord
    double precision :: detX, detY, detZ, dist

    logical :: ignite_detonation

!    integer :: point(MDIM)
!    double precision    :: dvol

    double precision :: dens, temp, eint, pres
    double precision :: xc12init, xne22init
    double precision :: flame, flamedot
    double precision :: phi_fa, phi_aq, phi_qn
    double precision :: ye, dyi_qn, dqbar_qn
    double precision :: qdot, edotnu
    double precision :: phi_fa_det

    double precision, parameter :: inflame_threshold = 1.e-6
    double precision :: cgsMeVperGram = 9.6485e17

    double precision    :: dti

    ! for saving some initial states
    double precision :: phi_fa_i, phi_aq_i, phi_qn_i
    double precision :: eint_i, pres_i
    double precision :: qbar_i, ye_i, yi_i, dye_n

    double precision :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a

    ! for calculating time evolution of properties
    double precision :: phi_aq_dot, phi_qn_dot
    double precision :: r_c
    double precision :: yi_q_est
    double precision :: tau_nsqe_inv, tau_nse_inv, maxExp, testExp  
    double precision :: expq

    double precision :: qbar

    ! for predictions of nse final state
    double precision :: s
    double precision :: qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d
    double precision :: qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p
    double precision :: qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot

    ! for calculation of reaction rate
    double precision :: nasigmav, t9, t9a

    ! for estimating reaction temperatures
    double precision :: qbar_unburned, yi_unburned, qbar_burned, yi_burned, x_bnf, qbar_noflame
    double precision :: yi_noflame, ye_noflame, eint_noflame

    type (eos_t) :: eos_state

100 format("*** igniting at ",f7.5,4(1x,es10.3)," ***")

!    if (autoDDT .and. thermalReact) then

       ! -------------------------------------------------------
       ! iterate over blocks and check for ignition conditions 
       ! -------------------------------------------------------
       ! determine the maximum number of ignition points to allow

!        ! initialize
!        ignition_num = 0
!        allocate(ignition_coords(dim*IgnNumMax),stat=istat)
!        if (istat/=0) call Driver_abortFlash("Unable to allocate ignition_coords")

!        ! store the number of ignition points in index zero
!        do n = 1, blockCount

!           spark_block = .false.
!           spark_block_radius = 0.e0
!           blockID = blockList(n)  

!           ! -------------------------------
!           ! 2.1 initialize quantities for this block
!           ! -------------------------------
!           call Grid_getBlkPtr(blockID,solnData)

!           ! get coordinate positions
!           call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!           xSizeCoord = blkLimitsGC(HIGH,IAXIS)
!           ySizeCoord = blkLimitsGC(HIGH,JAXIS)
!           zSizeCoord = blkLimitsGC(HIGH,KAXIS)
!           allocate(xCoord(xSizeCoord))
!           allocate(yCoord(ySizeCoord))
!           allocate(zCoord(zSizeCoord))
!           call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,xSizeCoord)
!           call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,ySizeCoord)
!           call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoord,zSizeCoord)

!           ! --------------------------------
!           ! 2.2 loop over all interior zones and check for ignition conditions 
!           ! --------------------------------
!           do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
!              do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
!                 do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

!                    ! --------------------------------
!                    ! 2.2.1 check for ignition conditions
!                    ! --------------------------------
! #ifndef FLAM_MSCALAR
!                    flame   = 0.e0
! #else
!                    flame   = solnData(FLAM_MSCALAR,i,j,k)
! #endif
!                    call paraSpark( xCoord(i), yCoord(j), zCoord(k), &
!                         solnData(DENS_VAR,i,j,k),        &
!                         solnData(PRES_VAR,i,j,k),        &
!                         solnData(PHFA_MSCALAR,i,j,k),    &
!                         flame,                           &
!                         solnData(CI_MSCALAR,i,j,k),      &
!                         solnData(NEI_MSCALAR,i,j,k),     &
!                         detX, detY, detZ, ignition_test )
!                    ! --------------------------------
!                    ! 2.2.2 store ignition coordinates
!                    ! --------------------------------
!                    if ( ignition_test ) then

!                       radius = sqrt( detX**2 + detY**2 + detZ**2 )

!                       ! check to see if we have already ignited on this block
!                       if (spark_block) then

!                          ! check which is further out from the center
!                          if (radius > spark_block_radius) then

!                             ! if current radius is larger, then overwrite previous
!                             ! ignition point
!                             spark_block_radius = radius
!                             ignition_coords((ignition_num-1)*dim + 1) = detX
! #if dim >= 2
!                             ignition_coords((ignition_num-1)*dim + 2) = detY
! #endif
! #if dim > 2
!                             ignition_coords((ignition_num-1)*dim + 3) = detZ
! #endif
!                          endif !! (radius > spark_block_radius)

!                          ! otherwise, our previous ignition point is further out
!                          ! and we will keep it

!                       else !! (spark_block)

!                          ! we have not already found an ignition condition on this
!                          ! block so, lets set a new one

!                          ! make sure we have enough space
!                          if ( ignition_num < IgnNumMax ) then

!                             spark_block = .true.
!                             ignition_coords( ignition_num*dim + 1 ) = detX
! #if dim >= 2
!                             ignition_coords( ignition_num*dim + 2 ) = detY
! #endif
! #if dim > 2
!                             ignition_coords( ignition_num*dim + 3 ) = detZ
! #endif
!                             ignition_num = ignition_num + 1

!                             ! we have too many ignition points, 
!                             ! need to create more space
!                          else !! ( ignition_num < IgnNumMax )

!                             call Driver_abortFlash("Not enough space to store all ignition points")

!                          endif !! ( ignition_num < IgnNumMax )

!                       endif !! (spark_block)

!                    endif !! (ignition_test)

!                 enddo
!              enddo
!           enddo

!           ! --------------------------------
!           ! 2.3 Finish up:
!           !     Release stuff
!           ! --------------------------------
!           call Grid_releaseBlkPtr(blockID,solnData)
!           deallocate(xCoord)
!           deallocate(yCoord)
!           deallocate(zCoord)
!        enddo

!        ! -------------------------------------------------------
!        ! 3.  communicate ignition points
!        ! -------------------------------------------------------
!        ! store number of ignition points in first index
!        if (ignition_num .eq. 0) then
!           deallocate(ignition_coords)
!        endif

!        ! allocate space for det coords
!        !     allocate(det_xCoord(IgnNumMax),stat=istat)
!        !     if (istat/=0) call Driver_abortFlash("Unable to allocate det_xCoord")
!        !     allocate(det_yCoord(IgnNumMax),stat=istat)
!        !     if (istat/=0) call Driver_abortFlash("Unable to allocate det_yCoord")
!        !     allocate(det_zCoord(IgnNumMax),stat=istat)
!        !     if (istat/=0) call Driver_abortFlash("Unable to allocate det_zCoord")

!        ! communicate detonation points and consolidate
!        ! this means only det points which are not near each other and
!        ! not near previous det points are returned, we still have to 
!        ! check these det points for being in ash
!        det_num = ignition_num
!        call paraAllSpark( ignition_coords, det_num, &
!             det_xCoord, det_yCoord, det_zCoord )

!        ! paraAllSpark has deallocated ignition_coords and has
!        ! allocated det_xCoord, det_yCoord, det_zCoord

!        ! -------------------------------------------------------
!        ! 4.  iterate over detonation points and blocks to
!        !     see if detonation points are valid
!        ! -------------------------------------------------------
!        ! make sure there are detonation points
!        if (det_num > 0) then

!           allocate(ignition_conditions(det_num),stat=istat)
!           if (istat/=0)  &
!                call Driver_abortFlash("Unable to allocate ignition_conditions")
!           ignition_conditions(:) = .true.

!           det_search:  do l = 1, det_num
!              ! ---------------------------
!              ! 4.1 iterate over all blocks
!              ! ---------------------------
!              do n = 1, blockCount

!                 blockID = blockList(n)

!                 ! -------------------------------
!                 ! 4.1.1 initialize quantities for this block
!                 ! -------------------------------
!                 call Grid_getBlkPtr(blockID,solnData)

!                 ! get coordinate positions, used for shock detection and
!                 ! non-flame burning proximity detection
!                 call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!                 xSizeCoord = blkLimitsGC(HIGH,IAXIS)
!                 ySizeCoord = blkLimitsGC(HIGH,JAXIS)
!                 zSizeCoord = blkLimitsGC(HIGH,KAXIS)
!                 allocate(xCoord(xSizeCoord))
!                 allocate(yCoord(ySizeCoord))
!                 allocate(zCoord(zSizeCoord))
!                 call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,xSizeCoord)
!                 call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,ySizeCoord)
!                 call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoord,zSizeCoord)

!                 ! --------------------------------
!                 ! 4.1.2 loop over all interior zones
!                 ! --------------------------------
!                 do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
!                    do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
!                       do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

!                          ! --------------------------------
!                          ! 4.1.2.1 check if we are near a detonation point
!                          ! --------------------------------
!                          dist = ( det_xCoord(l) - xCoord(i) )**2
! #if dim >= 2
!                          dist = dist + ( det_yCoord(l) - yCoord(j) )**2
! #endif
! #if dim > 2
!                          dist = dist + ( det_zCoord(l) - zCoord(k) )**2
! #endif
!                          dist = sqrt( dist )

!                          if (dist <= IgnRad) then

!                             ! check if we are in ash, and we are at ignition density
!                             if ( solnData(PHFA_MSCALAR,i,j,k) > IgnPhfa ) then

!                                ! we are in ash 
!                                ! so invalidate detonation point
!                                ignition_conditions(l) = .false.

!                                ! we are done so finish up with this block
!                                call Grid_releaseBlkPtr(blockID,solnData)
!                                deallocate(xCoord)
!                                deallocate(yCoord)
!                                deallocate(zCoord)

!                                ! move to the next detonation point
!                                cycle det_search
!                             endif
!                          endif

!                       enddo
!                    enddo
!                 enddo

!                 ! --------------------------------
!                 ! 4.1.3 Finish up:
!                 !     Release stuff
!                 ! --------------------------------
!                 call Grid_releaseBlkPtr(blockID,solnData)
!                 deallocate(xCoord)
!                 deallocate(yCoord)
!                 deallocate(zCoord)
!              enddo

!           enddo det_search

!           ! -----------------------
!           ! 4.2 communicate any invalidated detonation points
!           ! -----------------------
!           call paraAllIgnite( ignition_conditions, det_num )

!           ! announce detonation points and save to global list
!           if (IgnNum + count(ignition_conditions) < IgnNumMax) then
!              do l = 1, det_num
!                 if ( ignition_conditions(l) ) then
!                    if (meshMe .eq. MASTER_PE) then
!                       detX = det_xCoord(l)
! #if dim >= 2
!                       detY = det_yCoord(l)
! #else
!                       detY = 0.e0
! #endif
! #if dim > 2
!                       detZ = det_zCoord(l)
! #else
!                       detZ = 0.e0
! #endif
!                       open (burn_lun, file=detIgnFileName, &
!                            position='append', iostat=istat)
!                       if (istat/=0) call Driver_abortFlash("Unable to open detonation ignition points file")

!                       write (burn_lun,*) time, detX, detY, detZ

!                       close(unit=burn_lun)

!                    endif

!                    IgnNum = IgnNum + 1
!                    IgnTime( IgnNum ) = time
!                    IgnX( IgnNum ) = det_xCoord(l)
! #if dim >= 2
!                    IgnY( IgnNum ) = det_yCoord(l)
! #else
!                    IgnY( IgnNum ) = 0.e0
! #endif
! #if dim > 2
!                    IgnZ( IgnNum ) = det_zCoord(l)
! #else
!                    IgnZ( IgnNum ) = 0.e0
! #endif
!                 endif
!              enddo
!           else
!              call Driver_abortFlash("Not enough space to save detonation points")
!           endif

!        endif

!     endif

!       dx = xCoord(2)-xCoord(1) ! also assume square grid

       ! shock detect if burning is turned off in shocks
!       if (thermalReact .and. (.NOT. useShockBurn)) then
!          call Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, (/0,0,0/), &
!               xCoord,yCoord,zCoord)
!       else
!          shock(:,:,:) = 0
!       endif

!        !  Check for proximity of a reacting region for each cell
!        !  this is used to help control thermal burning inside flame
!        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
!           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
!              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
!                 !  proximity is in units of flame width (rounded up to nearest cell)
!                 ! no need to sqrt react_proximity because we are comparing to 1.0
!                 react_proximity(i,j,k) = 2.0  ! > 1 supresses reaction in flame
!                 ! allow burning with no flame module compiled in
!                 ! proximity doesn't matter with no flame, but be safe and just set leave it default
! #ifdef FLAM_MSCALAR
!                 if (thermalReact) then
!                    maxdi = int(ceiling(flamewidth/dx))
!                    do ii = -maxdi, maxdi
!                       do jj = -maxdi*K2D, maxdi*K2D
!                          do kk = -maxdi*K3D, maxdi*K3D
!                             if ( solnData(PHFA_MSCALAR,i+ii,j+jj,k+kk) - solnData(FLAM_MSCALAR,i+ii,j+jj,k+kk) &
!                                  > thermalReactInFlameThreshold ) then
!                                react_proximity(i,j,k) = min(react_proximity(i,j,k),dble(ii**2+jj**2+kk**2)/maxdi**2)
!                             endif
!                          enddo
!                       enddo
!                    enddo
!                 endif
! #endif
!              enddo
!           enddo
!        enddo

       ! --------------------------------
       ! set up internal energy, flame inputs
       ! and check detonation ignition points
       ! --------------------------------

       eint = state_in % e
       pres = state_in % p
       
       flame    = state_in % aux(1) ! UFLAM
       flamedot = state_in % aux(2) ! UFLDT
       
       ! ! --------------------------------
       ! ! check if we are near a detonation point
       ! ! --------------------------------

       ! ignite_detonation = .false.
       ! if (autoDDT) radius = IgnRad
       ! do l = 1, IgnNum

       !    dist = sqrt((IgnX(l) - xCoord(i))**2 + (IgnY(l) - yCoord(j))**2 + (IgnZ(l) - zCoord(k))**2)

       !    if (.not. autoDDT) radius = IgnR(l)
       !    if ( IgnTime(l) < time+dt .and. IgnTime(l) >= time &
       !         .and. dist <= radius ) then
       !       write (6,100) IgnTime(l), IgnX(l), IgnY(l), IgnZ(l), radius
       !       ignite_detonation = .true.
       !       phi_fa_det = 1.0e0
       !    endif
       ! enddo

       ! ! --------------------------------
       ! ! evolve progress variables and update NSE grid quantities
       ! ! --------------------------------

       ! if (shock == 1) then
       !    qdot = 0.0
       !    edotnu = 0.0
       !    return
       ! endif

       !--------------------------------------------------------
       ! initialize some information
       ! including getting unburned properties and saving some initial data
       !--------------------------------------------------------

       ! inverse timestep
       dti    = 1.e0/dt
       
       xc12init = state_in % xn(1)
       xne22init = state_in % xn(3)

       phi_fa = state_in % aux(4)
       phi_aq = state_in % aux(5)
       phi_qn = state_in % aux(6)
       dyi_qn = state_in % aux(8)       

       call paraFuelAshProperties(xc12init,xne22init, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

       ! truncate advection errors which might have pushed out of limits
       phi_fa = max(0.0,min(1.0,phi_fa))
       phi_aq = max(0.0,min(phi_fa,phi_aq))
       phi_qn = max(0.0,min(phi_aq,phi_qn))
       
       dyi_qn = max(0.0,min(1.0,dyi_qn))
       
       ! initial plasma properties
       qbar_i  = (1.0-phi_fa)*qbar_f + (phi_fa-phi_aq)*qbar_a + dqbar_qn
       ye_i  = ye
       yi_i  = (1.0-phi_fa)*yi_f + (phi_fa-phi_aq)*yi_a + dyi_qn
       
       dye_n = max( ye - (1.0-phi_qn)*ye_a, 0.0 )

       pres_i   = pres
       eint_i   = eint
       phi_fa_i = phi_fa
       phi_aq_i = phi_aq
       phi_qn_i = phi_qn

       !--------------------------------------------------------------------------
       ! Calculate properties (qbar, yi, temp) of NSE final state
       !--------------------------------------------------------------------------
       ! We predict the properties of the final NSE state from the current 
       ! (partially burned) state, assuming the burn will occur at either
       ! constant density (isochoric) or constant pressure (isobaric).
       ! The properties have been pre-calculated and tabulated for speed.  Inside a
       ! flame we use an interpolation between the isobaric and isochoric predictions
       ! In addtion to qbar, yi, and temp of the final state, we also obtain
       ! edot (neutrino losses) and yedot (neutronization rate) for that state.
       ! These are used in the evolution of partially burned zones.
       !--------------------------------------------------------------------------
       if ( flame >= 0.9999 .or. flame < inflame_threshold ) then
          ! isochoric prediction
          call NSE_finalAtDens(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
               ye_i, dens, eint_i - qbar_i*cgsMeVperGram)

       else if (flame >= 0.99) then
          ! crossover region, linearly average
          call NSE_finalAtDens(qbar_finnse_d, sumyi_finnse_d, temp_finnse_d, edot_d, yedot_d, &
               ye_i, dens, eint_i - qbar_i*cgsMeVperGram)
          call NSE_finalAtPres(qbar_finnse_p, sumyi_finnse_p, temp_finnse_p, edot_p, yedot_p, &
               ye_i, pres_i, eint_i + pres_i/dens - qbar_i*cgsMeVperGram)
          s = (flame - 0.99)/(0.9999-0.99)
          qbar_finnse = s*qbar_finnse_d + (1.0-s)*qbar_finnse_p
          sumyi_finnse = 1.0/(s/sumyi_finnse_d + (1.0-s)/sumyi_finnse_p)
          temp_finnse = s*temp_finnse_d + (1.0-s)*temp_finnse_p
          edot = s*edot_d + (1.0-s)*edot_p
          yedot = s*yedot_d + (1.0-s)*yedot_p

       else
          ! isobaric prediction
          call NSE_finalAtPres(qbar_finnse, sumyi_finnse, temp_finnse, edot, yedot, &
               ye_i, pres_i, eint_i + pres_i/dens - qbar_i*cgsMeVperGram)
       endif

       !-------------------------------------------------------------------
       !   calculate timescale calibrated in Calder etal 2007
       !-------------------------------------------------------------------

       testExp = 182.06e9/temp_finnse - 46.054e0
       ! protect from overflow for portability
       maxExp = log(HUGE(1.0))
       if (testExp .GE. maxExp) then
          tau_nsqe_inv = 0.0
       else
          tau_nsqe_inv = exp(-testExp)
       endif
       ! ajk     log(tau_nse)  = (dens**0.2e0)*exp(179.7e9/temp - 40.5e0)
       ! fang1   log(tau_nse)  = (dens**0.2e0)*exp(207.76e9/temp - 47.262e0)
       ! calder etal      log(tau_nse)  = exp(196.02e9/temp_finnse - 41.645e0)
       ! townsley etal '09 tau_nse = exp(-47.36 + 197.7/T_9)
       if ( flame >= 0.9999 .or. &
            ( thermalReact .and. &
            ( flame < inflame_threshold .or. (phi_fa > flame + 0.1) ) ) ) then
          ! use local temperature to get NSE timescale outside flame
          ! (particularly in detonations)
          !     testExp = 196.02e9/temp - 41.645e0
          !si     testExp = 197.7e9/temp - 47.36e0
          testExp = 201.0e9/temp - 46.77e0
       else
          ! in flame with no apparent thermal burning, use estimated final temperature
          ! this is just a rate, so a discontinuous crossover shouldn't be a big deal
          !     testExp = 196.02e9/temp_finnse - 41.645e0
          !si     testExp = 197.7e9/temp_finnse - 47.36e0
          testExp = 201.0e9/temp_finnse - 46.77e0
       endif
       if (testExp .GE. maxExp) then
          tau_nse_inv = 0.0
       else
          tau_nse_inv = exp(-testExp)
       endif


       !-----------------------------------------------------------------------
       !-----------------------------------------------------------------------

       ! if we are in the flame, only include thermal reactions if there is nearby
       ! thermal burning.  This prevents flame burning from causing thermal reactions
       ! without outside influence (such as an arriving detonation).
       if (thermalReact .and. &
            .not. ( flame >= inflame_threshold .and. react_proximity > 1.0 ) ) then
          ! away from flame or in flame near other reactions, use thermal reactions

          ! first need a temperature
          if (flame < inflame_threshold .or. phi_fa > (1.e0-inflame_threshold) ) then
             ! away from flame or with nearly fully burned material ...
             ! burn based on nuclear reaction rate determined by t9
             t9 = temp/1.e9
          else
             ! in flame but not yet fully burned...
             ! estimate temperature of material if flame had not burned

             ! catch up corner case in which the flame just started burning
             !   (prevents a divide by zero risk below)
             phi_fa = max(phi_fa,inflame_threshold)
             ! x_bnf = fraction of material which was burned, but not by the flame
             x_bnf = max(0.e0, phi_fa-flame)
             qbar_unburned =  qbar_f
             yi_unburned =  yi_f
             qbar_burned = ( qbar_i - qbar_unburned*(1.0-phi_fa) ) / phi_fa
             yi_burned = ( yi_i - yi_unburned*(1.0-phi_fa) ) / phi_fa
             qbar_noflame = qbar_unburned*(1.e0-x_bnf) + qbar_burned*x_bnf
             yi_noflame = yi_unburned*(1.e0-x_bnf) + yi_burned*x_bnf
             !ye_noflame   = 0.5*(1-x_bnf) + ye_oldash_i*x_bnf
             ! just use current ye instead
             !ye_noflame = ye_i
             ye_noflame = ye_f
             eint_noflame = eint_i - cgsMeVperGram*(qbar_i-qbar_noflame)

             if (yi_noflame < 0.e0) then
                write (6,*) ' neg sumy', phi_fa, flame, dqbar_qn
             endif
             ! find a floor energy for local material
             eos_state % rho  = dens
             eos_state % T    = small_temp
             eos_state % abar = 1.e0/yi_noflame
             eos_state % zbar = ye_noflame * eos_state % abar
             call eos(eos_input_rt, eos_state)
             if ( eint_noflame < eos_state % e ) then
                t9 = small_temp/1.e9
             else
                eos_state % e = eint_noflame
                eos_state % T = temp ! guess value
                call eos(eos_input_re, eos_state)
                t9 = eos_state % T / 1.e9
             endif
          endif

          ! thermally burn based on temperature just estimated
          ! From Caughlin & Fowler 1988, ADNDT 40 283
          ! note this is UNSCREENED, should be screened at some point
          t9a = t9/(1.0+0.0396*t9)
          nasigmav = 4.27e26*t9a**(5.0/6)*t9**(-1.5)*exp(-84.165*t9a**(-1.0/3.0)-2.12e-3*t9**3)
          ! integrate
          !     phi_fa_dot = dens*xc12init*(1.0-phi_fa)**2/12.0*nasigmav
          ! and add as evolution after flame
          r_c = dens*xc12init/12.0*nasigmav
          phi_fa = 1.0 - (1.0-phi_fa)/(1.0 + r_c*dt*(1.0-phi_fa))

          ! ignition override
          if (ignite_detonation) then
             ! fully burn this zone in this timestep (will hit limiter)
             phi_fa = phi_fa_det
          endif

       endif

       ! follow flame variable in any case
       ! (second operator split step if thermal burning is activated)
       phi_fa = max(0.0,min(1.e0,phi_fa+max(0.0e0,flamedot)*dt))


       !------------------------------------------------------
       ! 3 Update progress variables by direct integration, assuming
       !   phi_fa is constant (i.e. evolves much faster than others)

       expq = exp(-dt*tau_nsqe_inv)

       ! integration of     phi_aq_dot = (phi_fa-phi_aq)/tau_nsqe
       phi_aq = phi_fa - (phi_fa-phi_aq_i)*expq
       phi_aq = max(0.0,min(phi_fa,phi_aq))
       phi_aq_dot = (phi_aq-phi_aq_i)*dti

       ! integration of     phi_qn_dot = (phi_aq-phi_qn)^2/tau_nsqe
       ! assuming phi_aq is constant (i.e. operator split)
       ! need to handle case when phi_aq-phi_qn_i is zero
       if ((phi_aq-phi_qn_i) == 0.0) then
          phi_qn = phi_aq
       else
          phi_qn = phi_aq - 1.0/( 1.0/(phi_aq-phi_qn_i) + tau_nse_inv*dt )
       endif
       phi_qn = max(0.0,min(phi_aq,phi_qn))
       phi_qn_dot = (phi_qn-phi_qn_i)*dti

       ! now update partial properties
       ! operator split to mirror later treatment of dyi_qn
       dqbar_qn = dqbar_qn + phi_aq_dot*dt*236.537/28
       dqbar_qn = (phi_aq*236.537/28+phi_qn*(qbar_finnse-236.537/28))*(1.0-expq) + dqbar_qn*expq

       yi_q_est = yiion(iSi28)

       ! now update dyi_qn
       ! assuming phi_qn is constant (operator split again)
       ! two steps, first change due to change of phi_aq, then evolution
       !  better considering operator splitting of phi_qn and phi_aq evolution
       dyi_qn = dyi_qn + phi_aq_dot*dt*yi_q_est
       dyi_qn = (phi_aq*yi_q_est+phi_qn*(sumyi_finnse-yi_q_est))*(1.0-expq) + dyi_qn*expq

       ! implicit update not necessary, but use actual change in phi_qn from above
       !   ( this is like operator splitting it)
       dye_n = dye_n + ye_a*phi_qn_dot*dt + phi_qn*yedot*dt

       ! Construct final state
       qbar  = (1.0-phi_fa)*qbar_f + (phi_fa-phi_aq)*qbar_a + dqbar_qn
       ye    = (1.0-phi_qn)*ye_a + dye_n

       ! calculate energy release including neutrino losses
       ! and rest mass differences
       qdot = cgsMeVperGram*(qbar - qbar_i)*dti  &
            - phi_qn*( yedot*N_A*c_light*c_light*(m_p + m_e - m_n) + edot )

       edotnu = phi_qn*edot

       

       ! update energy
       
       state_out % e = state_in % e + qdot * dt

!       point(IAXIS) = i
!       point(JAXIS) = j
!       point(KAXIS) = k

!       neutLossThisProcStep = neutLossThisProcStep + state_out % rho * dvol * edotnu * dt

!     if (autoDDT .and. thermalReact) then

!        ! deallocate stuff
!        if (det_num > 0) then
!           deallocate(ignition_conditions)  
!           deallocate(det_xCoord)
! #if dim >= 2
!           deallocate(det_yCoord)
! #endif
! #if dim > 2
!           deallocate(det_zCoord)
! #endif
!        endif

!     endif

  end subroutine actual_burner  

end module actual_burner_module
