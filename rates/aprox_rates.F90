module aprox_rates_module

  ! These rate routines come from the public_aprox13/19/21.f90 files.
  ! Only those used in the networks we have in this repository are kept.
  ! We modify the calling sequence to take the tfactors as an argument.
  ! We comment out the density derivatives because we never evolve density
  ! in our reaction networks.

  use tfactors_module

  implicit none

  double precision, allocatable :: rv(:), tv(:), datn(:,:,:)
  double precision, allocatable :: rfdm(:),rfd0(:),rfd1(:),rfd2(:)
  double precision, allocatable :: tfdm(:),tfd0(:),tfd1(:),tfd2(:)

#ifdef CUDA
  attributes(managed) :: rv, tv, datn, rfdm, rfd0, rfd1, rfd2, tfdm, tfd0, tfd1, tfd2
#endif

  !$acc declare create(rv, tv, datn, rfdm, rfd0, rfd1, rfd2, tfdm, tfd0, tfd1, tfd2)

contains

  subroutine rates_init()

    implicit none

    integer :: j, k

    ! Allocate module arrays
    allocate(rv(6))
    allocate(tv(14))
    allocate(datn(2,6,14))
    allocate(rfdm(4))
    allocate(rfd0(4))
    allocate(rfd1(4))
    allocate(rfd2(4))
    allocate(tfdm(12))
    allocate(tfd0(12))
    allocate(tfd1(12))
    allocate(tfd2(12))

    rv = (/ 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 /)
    tv = (/ 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0 /)

    datn(1,:,:) = reshape( (/ -4.363, -3.091, -1.275, 1.073, 3.035, 4.825, &
			      -4.17, -2.964, -1.177, 1.085, 3.037, 4.826, &
			      -3.834, -2.727, -1.039, 1.104, 3.04, 4.826, &
			      -3.284, -2.418, -0.882, 1.129, 3.043, 4.827, &
			      -2.691, -2.093, -0.719, 1.159, 3.048, 4.827, &
			      -2.1675, -1.7668, -0.5573, 1.1947, 3.0527, 4.8272, &
			      -1.7095, -1.4462, -0.3991, 1.2358, 3.0577, 4.8276, &
        		      -1.3119, -1.1451, -0.2495, 1.2818, 3.0648, 4.8284, &
        		      -0.9812, -0.8612, -0.1084, 1.3336, 3.0738, 4.8295, &
        		      -0.682, -0.595, 0.028, 1.386, 3.084, 4.831, &
                              -0.4046, -0.3523, 0.1605, 1.4364, 3.0957, 4.8333, &
                              -0.1636, -0.1352, 0.2879, 1.4861, 3.1092, 4.8365, &
                              0.0461, 0.0595, 0.4105, 1.5354, 3.1242, 4.8405, &
                              0.2295, 0.235, 0.5289, 1.5842, 3.1405, 4.845 /), &
                           (/ 6, 14 /) )

    datn(2,:,:) = reshape( (/ -4.539, -3.097, -1.134, 1.525, 3.907, 6.078, &
        		      -4.199, -2.905, -1.024, 1.545, 3.91, 6.079, &
        		      -3.736, -2.602, -0.851, 1.578, 3.916, 6.08, &
        		      -3.052, -2.206, -0.636, 1.623, 3.923, 6.081, &
        		      -2.31, -1.766, -0.396, 1.678, 3.931, 6.082, &
        		      -1.6631, -1.319, -0.1438, 1.7471, 3.9409, 6.0829, &
        		      -1.1064, -0.8828, 0.1094, 1.8279, 3.9534, 6.0841, &
       			      -0.6344, -0.496, 0.3395, 1.9168, 3.9699, 6.0862, &
        		      -0.2568, -0.1555, 0.5489, 2.0163, 3.9906, 6.0893, &
        		      0.081, 0.158, 0.746, 2.114, 4.013, 6.093, &
        		      0.3961, 0.4448, 0.9304, 2.2026, 4.0363, 6.0976, &
        		      0.6673, 0.6964, 1.0985, 2.2849, 4.0614, 6.1033, &
        		      0.9009, 0.9175, 1.2525, 2.3619, 4.0882, 6.1099, &
        		      1.1032, 1.113, 1.3947, 2.4345, 4.1161, 6.1171 /), &
                           (/ 6, 14 /) )



    ! Evaluate the cubic interp parameters for ni56 electron capture
    ! which is used in the langanke subroutine.

    do k = 2, 4
       rfdm(k)=1./((rv(k-1)-rv(k))*(rv(k-1)-rv(k+1))*(rv(k-1)-rv(k+2)))
       rfd0(k)=1./((rv(k)-rv(k-1))*(rv(k)-rv(k+1))*(rv(k)-rv(k+2)))
       rfd1(k)=1./((rv(k+1)-rv(k-1))*(rv(k+1)-rv(k))*(rv(k+1)-rv(k+2)))
       rfd2(k)=1./((rv(k+2)-rv(k-1))*(rv(k+2)-rv(k))*(rv(k+2)-rv(k+1)))
    enddo

    do j = 2, 12
       tfdm(j)=1./((tv(j-1)-tv(j))*(tv(j-1)-tv(j+1))*(tv(j-1)-tv(j+2)))
       tfd0(j)=1./((tv(j)-tv(j-1))*(tv(j)-tv(j+1))*(tv(j)-tv(j+2)))
       tfd1(j)=1./((tv(j+1)-tv(j-1))*(tv(j+1)-tv(j))*(tv(j+1)-tv(j+2)))
       tfd2(j)=1./((tv(j+2)-tv(j-1))*(tv(j+2)-tv(j))*(tv(j+2)-tv(j+1)))
    enddo

  !$acc update device(rv, tv, datn, rfdm, rfd0, rfd1, rfd2, tfdm, tfd0, tfd1, tfd2)

  end subroutine rates_init



  subroutine rate_c12ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                        dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2, &
                        zz

    double precision, parameter :: q1 = 1.0d0/12.222016d0

    !$gpu

    ! c12(a,g)o16
    aa   = 1.0d0 + 0.0489d0*tf%t9i23
    daa  = -twoth*0.0489d0*tf%t9i53

    bb   = tf%t92*aa*aa
    dbb  = 2.0d0*(bb*tf%t9i + tf%t92*aa*daa)

    cc   = exp(-32.120d0*tf%t9i13 - tf%t92*q1)
    dcc  = cc * (oneth*32.120d0*tf%t9i43 - 2.0d0*tf%t9*q1)

    dd   = 1.0d0 + 0.2654d0*tf%t9i23
    ddd  = -twoth*0.2654d0*tf%t9i53

    ee   = tf%t92*dd*dd
    dee  = 2.0d0*(ee*tf%t9i + tf%t92*dd*ddd)

    ff   = exp(-32.120d0*tf%t9i13)
    dff  = ff * oneth*32.120d0*tf%t9i43

    gg   = 1.25d3 * tf%t9i32 * exp(-27.499*tf%t9i)
    dgg  = gg*(-1.5d0*tf%t9i + 27.499*tf%t9i2)

    hh   = 1.43d-2 * tf%t95 * exp(-15.541*tf%t9i)
    dhh  = hh*(5.0d0*tf%t9i + 15.541*tf%t9i2)

    zz   = 1.0d0/bb
    f1   = cc*zz
    df1  = (dcc - f1*dbb)*zz

    zz   = 1.0d0/ee
    f2   = ff*zz
    df2  = (dff - f2*dee)*zz

    term    = 1.04d8*f1  + 1.76d8*f2 + gg + hh
    dtermdt = 1.04d8*df1 + 1.76d8*df2 + dgg + dhh


    ! 1.7 times cf88 value
    term     = 1.7d0 * term
    dtermdt  = 1.7d0 * dtermdt

    fr    = term * den
    dfrdt = dtermdt * den * 1.0d-9
    !dfrdd = term

    rev    = 5.13d10 * tf%t932 * exp(-83.111*tf%t9i)
    drevdt = rev*(1.5d0*tf%t9i + 83.111*tf%t9i2)

    rr     = rev * term
    drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd  = 0.0d0

  end subroutine rate_c12ag

  ! This routine computes the nuclear reaction rate for 12C(a,g)16O and its inverse 
  ! using fit parameters from Deboer et al. 2017 (https://doi.org/10.1103/RevModPhys.89.035007).
  subroutine rate_c12ag_deboer17(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    type (tf_t)      :: tf

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd, &
			a0_nr,a1_nr,a2_nr,a3_nr,a4_nr,a5_nr,a6_nr, &
    			a0_r,a1_r,a2_r,a3_r,a4_r,a5_r,a6_r, &
                        term_a0_nr,term_a1_nr,term_a2_nr,term_a3_nr, &
                        term_a4_nr,term_a5_nr,term_a6_nr, &
                        term_a0_r,term_a1_r,term_a2_r,term_a3_r, &
                        term_a4_r,term_a5_r,term_a6_r, &
                        dterm_a0_nr,dterm_a1_nr,dterm_a2_nr,dterm_a3_nr,&
                        dterm_a4_nr,dterm_a5_nr,dterm_a6_nr, &
                        dterm_a0_r,dterm_a1_r,dterm_a2_r,dterm_a3_r,&
                        dterm_a4_r,dterm_a5_r,dterm_a6_r, &
                        term_nr,term_r,dterm_nr,dterm_r,  &
			term,dtermdt,rev,drevdt

    !$gpu
    
    ! from Table XXVI of deboer + 2017
    ! non-resonant contributions to the reaction
    a0_nr = 24.1d0
    a1_nr = 0d0 
    a2_nr = -32d0
    a3_nr = -5.9d0
    a4_nr = 1.8d0
    a5_nr = -0.17d0
    a6_nr = -twoth

    term_a0_nr = exp(a0_nr)
    term_a1_nr = exp(a1_nr*tf%t9i)
    term_a2_nr = exp(a2_nr*tf%t9i13)
    term_a3_nr = exp(a3_nr*tf%t913)
    term_a4_nr = exp(a4_nr*tf%t9)
    term_a5_nr = exp(a5_nr*tf%t953)
    term_a6_nr = tf%t9**a6_nr

    term_nr = term_a0_nr * term_a1_nr * term_a2_nr * &
              term_a3_nr * term_a4_nr * term_a5_nr * &
              term_a6_nr

    dterm_a0_nr = 0d0
    dterm_a1_nr = 0d0
    dterm_a2_nr = -a2_nr*tf%t9i43*term_a2_nr/3d0
    dterm_a3_nr = a3_nr*tf%t9i23*term_a3_nr/3d0
    dterm_a4_nr = a4_nr*term_a4_nr
    dterm_a5_nr = a5_nr*tf%t923*term_a5_nr*fiveth
    dterm_a6_nr = tf%t9i*a6_nr*tf%t9**a6_nr

    dterm_nr = (term_a0_nr * term_a1_nr * dterm_a2_nr * term_a3_nr * term_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * dterm_a3_nr * term_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * dterm_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * term_a4_nr * dterm_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * term_a4_nr * term_a5_nr * dterm_a6_nr)

    ! resonant contributions to the reaction
    a0_r = 7.4d0
    a1_r = -30d0
    a2_r = 0d0
    a3_r = 0d0
    a4_r = 0d0
    a5_r = 0d0
    a6_r = -3.0d0/2.0d0

    term_a0_r = exp(a0_r)
    term_a1_r = exp(a1_r*tf%t9i)
    term_a2_r = exp(a2_r*tf%t9i13)
    term_a3_r = exp(a3_r*tf%t913)
    term_a4_r = exp(a4_r*tf%t9)
    term_a5_r = exp(a5_r*tf%t953)
    term_a6_r = tf%t9**a6_r

    term_r = term_a0_r * term_a1_r * term_a2_r * &
              term_a3_r * term_a4_r * term_a5_r * &
              term_a6_r

    dterm_a0_r = 0d0
    dterm_a1_r = -a1_r*tf%t9i2*term_a1_r
    dterm_a2_r = 0d0
    dterm_a3_r = 0d0
    dterm_a4_r = 0d0
    dterm_a5_r = 0d0
    dterm_a6_r = tf%t9i*a6_r*tf%t9**a6_r
    
    dterm_r = (term_a0_r * dterm_a1_r * term_a6_r) + &
	      (term_a0_r * term_a1_r * dterm_a6_r)
    

    ! full rate is the sum of resonant and non-resonant contributions
    term = term_nr + term_r 
    dtermdt = dterm_nr + dterm_r

    fr    = term * den
    dfrdt = dtermdt * den * 1.0d-9

    ! first term is 9.8685d9 * T9**(2/3) * (M0*M1/M3)**(3/2) 
    ! see iliadis 2007 eqn. 3.44
    ! ratio of partition functions are assumed to be unity
    rev    = 5.1345573d10 * tf%t932 * exp(-83.114082*tf%t9i)
    drevdt = rev*(1.5d0*tf%t9i + 83.114082*tf%t9i2)

    rr     = rev * term
    drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9


  end subroutine rate_c12ag_deboer17


  subroutine rate_tripalf(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac, &
                        drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1

    double precision, parameter :: rc28   = 0.1d0
    double precision, parameter :: q1     = 1.0d0/0.009604d0
    double precision, parameter :: q2     = 1.0d0/0.055225d0

    !$gpu

    ! triple alfa to c12
    ! this is a(a,g)be8
    aa    = 7.40d+05 * tf%t9i32 * exp(-1.0663*tf%t9i)
    daa   = aa*(-1.5d0*tf%t9i  + 1.0663*tf%t9i2)

    bb    = 4.164d+09 * tf%t9i23 * exp(-13.49*tf%t9i13 - tf%t92*q1)
    dbb   = bb*(-twoth*tf%t9i + oneth*13.49*tf%t9i43 - 2.0d0*tf%t9*q1)

    cc    = 1.0d0 + 0.031*tf%t913 + 8.009*tf%t923 + 1.732*tf%t9 &
          + 49.883*tf%t943 + 27.426*tf%t953
    dcc   = oneth*0.031*tf%t9i23 + twoth*8.009*tf%t9i13 + 1.732 &
          + fourth*49.883*tf%t913 + fiveth*27.426*tf%t923

    r2abe    = aa + bb * cc
    dr2abedt = daa + dbb*cc + bb*dcc


    ! this is be8(a,g)c12
    dd    = 130.0d0 * tf%t9i32 * exp(-3.3364*tf%t9i)
    ddd   = dd*(-1.5d0*tf%t9i + 3.3364*tf%t9i2)

    ee    = 2.510d+07 * tf%t9i23 * exp(-23.57*tf%t9i13 - tf%t92*q2)
    dee   = ee*(-twoth*tf%t9i + oneth*23.57*tf%t9i43 - 2.0d0*tf%t9*q2)

    ff    = 1.0d0 + 0.018*tf%t913 + 5.249*tf%t923 + 0.650*tf%t9 + &
         19.176*tf%t943 + 6.034*tf%t953
    dff   = oneth*0.018*tf%t9i23 + twoth*5.249*tf%t9i13 + 0.650 &
          + fourth*19.176*tf%t913 + fiveth*6.034*tf%t923

    rbeac    = dd + ee * ff
    drbeacdt = ddd + dee * ff + ee * dff


    ! a factor
    xx    = rc28 * 1.35d-07 * tf%t9i32 * exp(-24.811*tf%t9i)
    dxx   = xx*(-1.5d0*tf%t9i + 24.811*tf%t9i2)


    ! high temperature rate
    if (tf%t9.gt.0.08) then
       term    = 2.90d-16 * r2abe * rbeac + xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac &
                 + 2.90d-16 * r2abe * drbeacdt &
                 + dxx


    ! low temperature rate
    else
       uu   = 0.8d0*exp(-(0.025*tf%t9i)**3.263)
       yy   = 0.2d0 + uu
       ! fxt yy   = 0.01 + 0.2d0 + uu
       dyy  = uu * 3.263*(0.025*tf%t9i)**2.263 * (0.025*tf%t9i2)
       vv   = 4.0d0*exp(-(tf%t9/0.025)**9.227)
       zz   = 1.0d0 + vv
       dzz  = vv * 9.227*(tf%t9/0.025)**8.227 * 40.0d0
       aa   = 1.0d0/zz
       f1   = 0.01d0 + yy * aa
       ! fxt f1   = yy * aa
       df1  = (dyy - f1*dzz)*aa
       term = 2.90d-16 * r2abe * rbeac * f1 +  xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac * f1 &
                 + 2.90d-16 * r2abe * drbeacdt * f1 &
                 + 2.90d-16 * r2abe * rbeac * df1 &
                 + dxx
    end if


    ! rates
    !      term    = 1.2d0 * term
    !      dtermdt = 1.2d0 * term

    fr    = term * den * den
    dfrdt = dtermdt * den * den * 1.0d-9
    !dfrdd = 2.0d0 * term * den

    rev    = 2.00d+20*tf%t93*exp(-84.424*tf%t9i)
    drevdt = rev*(3.0d0*tf%t9i + 84.424*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_tripalf



  subroutine rate_c12c12(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                        aa,zz

    !$gpu

    ! c12 + c12 reaction
    aa      = 1.0d0 + 0.0396*tf%t9
    zz      = 1.0d0/aa

    t9a     = tf%t9*zz
    dt9a    = (1.0d0 -  t9a*0.0396)*zz

    zz      = dt9a/t9a
    t9a13   = t9a**oneth
    dt9a13  = oneth*t9a13*zz

    t9a56   = t9a**fivsix
    dt9a56  = fivsix*t9a56*zz

    term    = 4.27d+26 * t9a56 * tf%t9i32 * &
         exp(-84.165/t9a13 - 2.12d-03*tf%t93)
    dtermdt = term*(dt9a56/t9a56 - 1.5d0*tf%t9i &
            + 84.165/t9a13**2*dt9a13 - 6.36d-3*tf%t92)

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rr    = 0.0d0
    drrdt = 0.0d0
    !drrdd = 0.0d0

  end subroutine rate_c12c12



  subroutine rate_c12o16(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23, &
                        t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz

    !$gpu

    ! c12 + o16 reaction; see cf88 references 47-4
    if (tf%t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*tf%t9
       zz     = 1.0d0/aa

       t9a    = tf%t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa      = exp(-0.18*t9a*t9a)
       daa     = -aa * 0.36 * t9a * dt9a

       bb      = 1.06d-03*exp(2.562*t9a23)
       dbb     = bb * 2.562 * dt9a23

       cc      = aa + bb
       dcc     = daa + dbb

       zz      = 1.0d0/cc
       term    = 1.72d+31 * t9a56 * tf%t9i32 * exp(-106.594/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5d0*tf%t9i &
                       + 106.594/t9a23*dt9a13 - zz*dcc)

    else
       ! term    = 2.6288035d-29
       term    = 0.0d0
       dtermdt = 0.0d0
    endif


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rr    = 0.0d0
    drrdt = 0.0d0
    !drrdd = 0.0d0

  end subroutine rate_c12o16



  subroutine rate_o16o16(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt

    !$gpu

    ! o16 + o16
    term  = 7.10d36 * tf%t9i23 * &
         exp(-135.93 * tf%t9i13 - 0.629*tf%t923 &
         - 0.445*tf%t943 + 0.0103*tf%t9*tf%t9)

    dtermdt = -twoth*term*tf%t9i &
         + term * (oneth*135.93*tf%t9i43 - twoth*0.629*tf%t9i13 &
         - fourth*0.445*tf%t913 + 0.0206*tf%t9)


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rr    = 0.0d0
    drrdt = 0.0d0
    !drrdd = 0.0d0

  end subroutine rate_o16o16



  subroutine rate_o16ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                        cc,dcc,term2,dterm2,rev,drevdt

    double precision, parameter :: q1 = 1.0d0/2.515396d0

    !$gpu

    ! o16(a,g)ne20
    term1   = 9.37d9 * tf%t9i23 * exp(-39.757*tf%t9i13 - tf%t92*q1)
    dterm1  = term1*(-twoth*tf%t9i + oneth*39.757*tf%t9i43 - 2.0d0*tf%t9*q1)

    aa      = 62.1 * tf%t9i32 * exp(-10.297*tf%t9i)
    daa     = aa*(-1.5d0*tf%t9i + 10.297*tf%t9i2)

    bb      = 538.0d0 * tf%t9i32 * exp(-12.226*tf%t9i)
    dbb     = bb*(-1.5d0*tf%t9i + 12.226*tf%t9i2)

    cc      = 13.0d0 * tf%t92 * exp(-20.093*tf%t9i)
    dcc     = cc*(2.0d0*tf%t9i + 20.093*tf%t9i2)

    term2   = aa + bb + cc
    dterm2  = daa + dbb + dcc

    term    = term1 + term2
    dtermdt = dterm1 + dterm2


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 5.65d+10*tf%t932*exp(-54.937*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 54.937*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_o16ag



  subroutine rate_ne20ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                        term2,dterm2,term3,dterm3,rev,drevdt,zz

    double precision, parameter :: rc102 = 0.1d0
    double precision, parameter :: q1    = 1.0d0/4.923961d0

    !$gpu

    ! ne20(a,g)mg24
    aa   = 4.11d+11 * tf%t9i23 * exp(-46.766*tf%t9i13 - tf%t92*q1)
    daa  = aa*(-twoth*tf%t9i + oneth*46.766*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb   = 1.0d0 + 0.009*tf%t913 + 0.882*tf%t923 + 0.055*tf%t9 &
         + 0.749*tf%t943 + 0.119*tf%t953
    dbb  = oneth*0.009*tf%t9i23 + twoth*0.882*tf%t9i13 + 0.055 &
         + fourth*0.749*tf%t913 + fiveth*0.119*tf%t923

    term1  = aa * bb
    dterm1 = daa * bb + aa * dbb


    aa   = 5.27d+03 * tf%t9i32 * exp(-15.869*tf%t9i)
    daa  = aa*(-1.5d0*tf%t9i + 15.869*tf%t9i2)

    bb   = 6.51d+03 * tf%t912 * exp(-16.223*tf%t9i)
    dbb  = bb*(0.5d0*tf%t9i + 16.223*tf%t9i2)

    term2  = aa + bb
    dterm2 = daa + dbb


    aa   = 42.1 * tf%t9i32 * exp(-9.115*tf%t9i)
    daa  = aa*(-1.5d0*tf%t9i + 9.115*tf%t9i2)

    bb   =  32.0 * tf%t9i23 * exp(-9.383*tf%t9i)
    dbb  = bb*(-twoth*tf%t9i + 9.383*tf%t9i2)

    term3  = rc102 * (aa + bb)
    dterm3 = rc102 * (daa + dbb)


    aa  = 5.0d0*exp(-18.960*tf%t9i)
    daa = aa*18.960*tf%t9i2

    bb  = 1.0d0 + aa
    dbb = daa

    zz      = 1.0d0/bb
    term    = (term1 + term2 + term3)*zz
    dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 6.01d+10 * tf%t932 * exp(-108.059*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 108.059*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_ne20ag



  subroutine rate_mg24ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,gg,dgg,hh,hhi,rev,drevdt

    double precision, parameter :: rc121 = 0.1d0

    !$gpu

    ! 24mg(a,g)28si
    aa    = 4.78d+01 * tf%t9i32 * exp(-13.506*tf%t9i)
    daa   = aa*(-1.5d0*tf%t9i + 13.506*tf%t9i2)

    bb    =  2.38d+03 * tf%t9i32 * exp(-15.218*tf%t9i)
    dbb   = bb*(-1.5d0*tf%t9i + 15.218*tf%t9i2)

    cc    = 2.47d+02 * tf%t932 * exp(-15.147*tf%t9i)
    dcc   = cc*(1.5d0*tf%t9i + 15.147*tf%t9i2)

    dd    = rc121 * 1.72d-09 * tf%t9i32 * exp(-5.028*tf%t9i)
    ddd   = dd*(-1.5d0*tf%t9i + 5.028*tf%t9i2)

    ee    = rc121* 1.25d-03 * tf%t9i32 * exp(-7.929*tf%t9i)
    dee   = ee*(-1.5d0*tf%t9i + 7.929*tf%t9i2)

    ff    = rc121 * 2.43d+01 * tf%t9i * exp(-11.523*tf%t9i)
    dff   = ff*(-tf%t9i + 11.523*tf%t9i2)

    gg    = 5.0d0*exp(-15.882*tf%t9i)
    dgg   = gg*15.882*tf%t9i2

    hh    = 1.0d0 + gg
    hhi   = 1.0d0/hh

    term    = (aa + bb + cc + dd + ee + ff) * hhi
    dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 6.27d+10 * tf%t932 * exp(-115.862*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 115.862*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_mg24ag



  subroutine rate_mg24ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,gg,dgg,term1,dterm1,term2,dterm2, &
                        rev,drevdt

    double precision, parameter :: rc148 = 0.1d0
    double precision, parameter :: q1    = 1.0d0/0.024649d0

    !$gpu

    ! 24mg(a,p)al27
    aa     = 1.10d+08 * tf%t9i23 * exp(-23.261*tf%t9i13 - tf%t92*q1)
    daa    = -twoth*aa*tf%t9i + aa*(23.261*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb     =  1.0d0 + 0.018*tf%t913 + 12.85*tf%t923 + 1.61*tf%t9 &
         + 89.87*tf%t943 + 28.66*tf%t953
    dbb    = oneth*0.018*tf%t9i23 + twoth*12.85*tf%t9i13 + 1.61 &
           + fourth*89.87*tf%t913 + fiveth*28.66*tf%t923

    term1  = aa * bb
    dterm1 = daa * bb + aa * dbb

    aa     = 129.0d0 * tf%t9i32 * exp(-2.517*tf%t9i)
    daa    = -1.5d0*aa*tf%t9i + aa*2.517*tf%t9i2

    bb     = 5660.0d0 * tf%t972 * exp(-3.421*tf%t9i)
    dbb    = 3.5d0*bb*tf%t9i +  bb*3.421*tf%t9i2

    cc     = rc148 * 3.89d-08 * tf%t9i32 * exp(-0.853*tf%t9i)
    dcc    = -1.5d0*cc*tf%t9i + cc*0.853*tf%t9i2

    dd     = rc148 * 8.18d-09 * tf%t9i32 * exp(-1.001*tf%t9i)
    ddd    = -1.5d0*dd*tf%t9i + dd*1.001*tf%t9i2

    term2  = aa + bb + cc + dd
    dterm2 = daa + dbb + dcc + ddd

    ee     = oneth*exp(-9.792*tf%t9i)
    dee    = ee*9.792*tf%t9i2

    ff     =  twoth * exp(-11.773*tf%t9i)
    dff    = ff*11.773*tf%t9i2

    gg     = 1.0d0 + ee + ff
    dgg    = dee + dff

    term    = (term1 + term2)/gg
    dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


    ! the rates
    rev      = 1.81 * exp(-18.572*tf%t9i)
    drevdt   = rev*18.572*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    !drrdd = term

  end subroutine rate_mg24ap



  subroutine rate_al27pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                        dd,ddd,ee,dee,ff,dff,gg,dgg

    !$gpu

    ! al27(p,g)si28
    ! champagne 1996

    aa  = 1.32d+09 * tf%t9i23 * exp(-23.26*tf%t9i13)
    daa = aa*(-twoth*tf%t9i + oneth*23.26*tf%t9i43)

    bb  = 3.22d-10 * tf%t9i32 * exp(-0.836*tf%t9i)*0.17
    dbb = bb*(-1.5d0*tf%t9i + 0.836*tf%t9i2)

    cc  = 1.74d+00 * tf%t9i32 * exp(-2.269*tf%t9i)
    dcc = cc*(-1.5d0*tf%t9i + 2.269*tf%t9i2)

    dd  = 9.92d+00 * tf%t9i32 * exp(-2.492*tf%t9i)
    ddd = dd*(-1.5d0*tf%t9i + 2.492*tf%t9i2)

    ee  = 4.29d+01 * tf%t9i32 * exp(-3.273*tf%t9i)
    dee = ee*(-1.5d0*tf%t9i + 3.273*tf%t9i2)

    ff  = 1.34d+02 * tf%t9i32 * exp(-3.654*tf%t9i)
    dff = ff*(-1.5d0*tf%t9i + 3.654*tf%t9i2)

    gg  = 1.77d+04 * (tf%t9**0.53) * exp(-4.588*tf%t9i)
    dgg = gg*(0.53*tf%t9i + 4.588*tf%t9i2)

    term    = aa + bb + cc + dd + ee + ff + gg
    dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.13d+11 * tf%t932 * exp(-134.434*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 134.434*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_al27pg



  subroutine rate_al27pg_old(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,gg,dgg,hh,dhh,xx,dxx,yy,dyy,zz,dzz,pp, &
                        rev,drevdt

    double precision, parameter :: rc147 = 0.1d0
    double precision, parameter :: q1    = 1.0d0/0.024025d0

    !$gpu

    ! 27al(p,g)si28  cf88
    aa  = 1.67d+08 * tf%t9i23 * exp(-23.261*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*23.261*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb  = 1.0d0 + 0.018*tf%t913 + 5.81*tf%t923 + 0.728*tf%t9 &
         + 27.31*tf%t943 + 8.71*tf%t953
    dbb = oneth*0.018*tf%t9i23 + twoth*5.81*tf%t9i13 + 0.728 &
         + fourth*27.31*tf%t913 + fiveth*8.71*tf%t923

    cc  = aa*bb
    dcc = daa*bb + aa*dbb

    dd  = 2.20d+00 * tf%t9i32 * exp(-2.269*tf%t9i)
    ddd = dd*(-1.5d0*tf%t9i + 2.269*tf%t9i2)

    ee  = 1.22d+01 * tf%t9i32 * exp(-2.491*tf%t9i)
    dee = ee*(-1.5d0*tf%t9i + 2.491*tf%t9i2)

    ff  =  1.50d+04 * tf%t9 * exp(-4.112*tf%t9i)
    dff = ff*(tf%t9i + 4.112*tf%t9i2)

    gg  = rc147 * 6.50d-10 * tf%t9i32 * exp(-0.853*tf%t9i)
    dgg = gg*(-1.5d0*tf%t9i + 0.853*tf%t9i2)

    hh  = rc147 * 1.63d-10 * tf%t9i32 * exp(-1.001*tf%t9i)
    dhh = hh*(-1.5d0*tf%t9i + 1.001*tf%t9i2)

    xx     = oneth*exp(-9.792*tf%t9i)
    dxx    = xx*9.792*tf%t9i2

    yy     =  twoth * exp(-11.773*tf%t9i)
    dyy    = yy*11.773*tf%t9i2

    zz     = 1.0d0 + xx + yy
    dzz    = dxx + dyy

    pp      = 1.0d0/zz
    term    = (cc + dd + ee + ff + gg + hh)*pp
    dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.13d+11*tf%t932*exp(-134.434*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 134.434*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_al27pg_old



  subroutine rate_si28ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! si28(a,g)s32
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
    if (z .eq. 10.0) then
       daa = 0
    else
       daa   = 6.340d-2 + 2.0d0*2.541d-3*tf%t9 - 3.0d0*2.900d-4*tf%t92
    end if

    term    = 4.82d+22 * tf%t9i23 * exp(-61.015 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 61.015*tf%t9i13*(oneth*tf%t9i*aa - daa))
  
    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 6.461d+10 * tf%t932 * exp(-80.643*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 80.643*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_si28ag



  subroutine rate_si28ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! si28(a,p)p31
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 2.798d-3 + 2.0d0*2.763d-3*tf%t9 - 3.0d0*2.341d-4*tf%t92
    end if

    term    = 4.16d+13 * tf%t9i23 * exp(-25.631 * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*25.631*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 0.5825d0 * exp(-22.224*tf%t9i)
    drevdt   = rev*22.224*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    !drrdd = term

  end subroutine rate_si28ap



  subroutine rate_p31pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! p31(p,g)s32
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.928d-1 - 2.0d0*1.540d-2*tf%t9 + 3.0d0*6.444d-4*tf%t92
    end if

    term    = 1.08d+16 * tf%t9i23 * exp(-27.042 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 27.042*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 3.764d+10 * tf%t932 * exp(-102.865*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 102.865*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_p31pg



  subroutine rate_s32ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! s32(a,g)ar36
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 4.913d-2 + 2.0d0*4.637d-3*tf%t9 - 3.0d0*4.067d-4*tf%t92
    end if

    term    = 1.16d+24 * tf%t9i23 * exp(-66.690 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 66.690*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 6.616d+10 * tf%t932 * exp(-77.080*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 77.080*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_s32ag



  subroutine rate_s32ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! s32(a,p)cl35
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
    if (z .eq. 10) then
       daa = 0.0d0
    else
       daa   = 1.041d-1 - 2.0d0*1.368d-2*tf%t9 + 3.0d0*6.969d-4*tf%t92
    end if

    term    = 1.27d+16 * tf%t9i23 * exp(-31.044 * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*31.044*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 1.144 * exp(-21.643*tf%t9i)
    drevdt   = rev*21.643*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    !drrdd = term

  end subroutine rate_s32ap



  subroutine rate_cl35pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt

    !$gpu

    ! cl35(p,g)ar36
    aa    = 1.0d0 + 1.761d-1*tf%t9 - 1.322d-2*tf%t92 + 5.245d-4*tf%t93
    daa   = 1.761d-1 - 2.0d0*1.322d-2*tf%t9 + 3.0d0*5.245d-4*tf%t92
  

    term    =  4.48d+16 * tf%t9i23 * exp(-29.483 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 29.483*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 7.568d+10*tf%t932*exp(-98.722*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 98.722*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_cl35pg


  subroutine rate_ar36ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ar36(a,g)ca40
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.458d-1 - 2.0d0*1.069d-2*tf%t9 + 3.0d0*3.790d-4*tf%t92
    end if

    term    = 2.81d+30 * tf%t9i23 * exp(-78.271 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 78.271*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 6.740d+10 * tf%t932 * exp(-81.711*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 81.711*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_ar36ag



  subroutine rate_ar36ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ar36(a,p)k39
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 4.826d-3 - 2.0d0*5.534d-3*tf%t9 + 3.0d0*4.021d-4*tf%t92
    end if

    term    = 2.76d+13 * tf%t9i23 * exp(-34.922 * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*34.922*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 1.128*exp(-14.959*tf%t9i)
    drevdt   = rev*14.959*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    !drrdd = term

  end subroutine rate_ar36ap



  subroutine rate_k39pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! k39(p,g)ca40
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
    if (z .eq. 10) then
       daa = 0.0d0
    else
       daa   = 1.622d-1 - 2.0d0*1.119d-2*tf%t9 + 3.0d0*3.910d-4*tf%t92
    end if

    term    = 4.09d+16 * tf%t9i23 * exp(-31.727 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 31.727*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 7.600d+10 * tf%t932 * exp(-96.657*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 96.657*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_k39pg



  subroutine rate_ca40ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ca40(a,g)ti44
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.650d-2 + 2.0d0*5.973d-3*tf%t9 - 3.0d0*3.889d-4*tf%t92
    end if

    term    = 4.66d+24 * tf%t9i23 * exp(-76.435 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 76.435*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 6.843d+10 * tf%t932 * exp(-59.510*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 59.510*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_ca40ag



  subroutine rate_ca40ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ca40(a,p)sc43
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = -1.206d-2 + 2.0d0*7.753d-3*tf%t9 - 3.0d0*5.071d-4*tf%t92
    end if

    term    = 4.54d+14 * tf%t9i23 * exp(-32.177 * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*32.177*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 2.229 * exp(-40.966*tf%t9i)
    drevdt   = rev*40.966*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    !drrdd = term

  end subroutine rate_ca40ap



  subroutine rate_sc43pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! sc43(p,g)ca40
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.023d-1 - 2.0d0*2.242d-3*tf%t9 - 3.0d0*5.463d-5*tf%t92
    end if

    term    = 3.85d+16 * tf%t9i23 * exp(-33.234 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 33.234*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.525d+11 * tf%t932 * exp(-100.475*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 100.475*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_sc43pg



  subroutine rate_ti44ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ti44(a,g)cr48
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.066d-1 - 2.0d0*1.102d-2*tf%t9 + 3.0d0*5.324d-4*tf%t92
    end if

    term    = 1.37d+26 * tf%t9i23 * exp(-81.227 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 81.227*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 6.928d+10*tf%t932*exp(-89.289*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 89.289*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_ti44ag



  subroutine rate_ti44ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ti44(a,p)v47
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 2.655d-2 - 2.0d0*3.947d-3*tf%t9 + 3.0d0*2.522d-4*tf%t92
    end if

    term    = 6.54d+20 * tf%t9i23 * exp(-66.678 * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*66.678*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 1.104 * exp(-4.723*tf%t9i)
    drevdt   = rev*4.723*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    !drrdd = term

  end subroutine rate_ti44ap



  subroutine rate_v47pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! v47(p,g)cr48
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 9.979d-2 - 2.0d0*2.269d-3*tf%t9 - 3.0d0*6.662d-5*tf%t92
    end if

    term    = 2.05d+17 * tf%t9i23 * exp(-35.568 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 35.568*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 7.649d+10*tf%t932*exp(-93.999*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 93.999*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_v47pg



  subroutine rate_cr48ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! cr48(a,g)fe52
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 6.325d-2 - 2.0d0*5.671d-3*tf%t9 + 3.0d0*2.848d-4*tf%t92
    end if

    term    = 1.04d+23 * tf%t9i23 * exp(-81.420 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 81.420*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 7.001d+10 * tf%t932 * exp(-92.177*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 92.177*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_cr48ag



  subroutine rate_cr48ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! cr48(a,p)mn51
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.384d-2 + 2.0d0*1.081d-3*tf%t9 - 3.0d0*5.933d-5*tf%t92
    end if

    term    = 1.83d+26 * tf%t9i23 * exp(-86.741 * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*86.741*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 0.6087*exp(-6.510*tf%t9i)
    drevdt   = rev*6.510*tf%t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = rev * term

  end subroutine rate_cr48ap



  subroutine rate_mn51pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! mn51(p,g)fe52
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 8.922d-2 - 2.0d0*1.256d-3*tf%t9 - 3.0d0*9.453d-5*tf%t92
    end if

    term    = 3.77d+17 * tf%t9i23 * exp(-37.516 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 37.516*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.150d+11*tf%t932*exp(-85.667*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 85.667*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_mn51pg



  subroutine rate_fe52ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! fe52(a,g)ni56
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 7.846d-2 - 2.0d0*7.430d-3*tf%t9 + 3.0d0*3.723d-4*tf%t92
    end if

    term    = 1.05d+27 * tf%t9i23 * exp(-91.674 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 91.674*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 7.064d+10*tf%t932*exp(-92.850*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 92.850*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_fe52ag



  subroutine rate_fe52ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! fe52(a,p)co55
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.367d-2 + 2.0d0*7.428d-4*tf%t9 - 3.0d0*3.050d-5*tf%t92
    end if

    term    = 1.30d+27 * tf%t9i23 * exp(-91.674 * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*91.674*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 0.4597*exp(-9.470*tf%t9i)
    drevdt   = rev*9.470*tf%t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = rev * term

  end subroutine rate_fe52ap



  subroutine rate_co55pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! co55(p,g)ni56
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 9.894d-2 - 2.0d0*3.131d-3*tf%t9 - 3.0d0*2.160d-5*tf%t92
    end if

    term    = 1.21d+18 * tf%t9i23 * exp(-39.604 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 39.604*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.537d+11*tf%t932*exp(-83.382*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 83.382*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_co55pg



  subroutine rate_pp(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,aa,daa,bb,dbb

    !$gpu

    ! p(p,e+nu)d
    if (tf%t9 .le. 3.0) then
       aa   = 4.01d-15 * tf%t9i23 * exp(-3.380d0*tf%t9i13)
       daa  = aa*(-twoth*tf%t9i + oneth*3.380d0*tf%t9i43)

       bb   = 1.0d0 + 0.123d0*tf%t913 + 1.09d0*tf%t923 + 0.938d0*tf%t9
       dbb  = oneth*0.123d0*tf%t9i23 + twoth*1.09d0*tf%t9i13 + 0.938d0

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

    else
       term    = 1.1581136d-15
       dtermdt = 0.0d0
    end if

    ! rate
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rr    = 0.0d0
    drrdt = 0.0d0
    !drrdd = 0.0d0

  end subroutine rate_pp



  subroutine rate_png(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa

    !$gpu

    ! p(n,g)d
    ! smith,kawano,malany 1992

    aa      = 1.0d0 - 0.8504*tf%t912 + 0.4895*tf%t9 &
         - 0.09623*tf%t932 + 8.471e-3*tf%t92 &
         - 2.80e-4*tf%t952

    daa     =  -0.5d0*0.8504*tf%t9i12 + 0.4895 &
         - 1.5d0*0.09623*tf%t912 + 2.0d0*8.471e-3*tf%t9 &
         - 2.5d0*2.80e-4*tf%t932

    term    = 4.742e4 * aa
    dtermdt = 4.742e4 * daa


    ! wagoner,schramm 1977
    !      aa      = 1.0d0 - 0.86*tf%t912 + 0.429*tf%t9
    !      daa     =  -0.5d0*0.86*tf%t9i12 + 0.429

    !      term    = 4.4d4 * aa
    !      dtermdt = 4.4d4 * daa



    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 4.71d+09 * tf%t932 * exp(-25.82*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 25.82*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_png



  subroutine rate_dpg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb

    !$gpu

    ! d(p,g)he3
    aa      = 2.24d+03 * tf%t9i23 * exp(-3.720*tf%t9i13)
    daa     = aa*(-twoth*tf%t9i + oneth*3.720*tf%t9i43)

    bb      = 1.0d0 + 0.112*tf%t913 + 3.38*tf%t923 + 2.65*tf%t9
    dbb     = oneth*0.112*tf%t9i23 + twoth*3.38*tf%t9i13 + 2.65

    term    = aa * bb
    dtermdt = daa * bb + aa * dbb


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.63d+10 * tf%t932 * exp(-63.750*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 63.750*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_dpg



  subroutine rate_he3ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt

    !$gpu

    ! he3(n,g)he4
    term    = 6.62 * (1.0d0 + 905.0*tf%t9)
    dtermdt = 5.9911d3

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 2.61d+10 * tf%t932 * exp(-238.81*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 238.81*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_he3ng



  subroutine rate_he3he3(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb

    !$gpu

    ! he3(he3,2p)he4
    aa   = 6.04d+10 * tf%t9i23 * exp(-12.276*tf%t9i13)
    daa  = aa*(-twoth*tf%t9i + oneth*12.276*tf%t9i43)

    bb   = 1.0d0 + 0.034*tf%t913 - 0.522*tf%t923 - 0.124*tf%t9 &
         + 0.353*tf%t943 + 0.213*tf%t953
    dbb  = oneth*0.034*tf%t9i23 - twoth*0.522*tf%t9i13 - 0.124 &
         + fourth*0.353*tf%t913 + fiveth*0.213*tf%t923

    term    = aa * bb
    dtermdt = daa*bb + aa*dbb

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 3.39e-10 * tf%t9i32 * exp(-149.230*tf%t9i)
    drevdt   = rev*(-1.5d0*tf%t9i + 149.230*tf%t9i2)

    rr    = den * den * rev * term
    drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 2.0d0 * den * rev * term

  end subroutine rate_he3he3



  subroutine rate_he3he4(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,t9a,dt9a, &
                        t9a13,dt9a13,t9a56,dt9a56,zz

    !$gpu

    ! he3(he4,g)be7
    aa      = 1.0d0 + 0.0495*tf%t9
    daa     = 0.0495

    zz      = 1.0d0/aa
    t9a     = tf%t9*zz
    dt9a    = (1.0d0 - t9a*daa)*zz

    zz      = dt9a/t9a
    t9a13   = t9a**oneth
    dt9a13  = oneth*t9a13*zz

    t9a56   = t9a**fivsix
    dt9a56  = fivsix*t9a56*zz

    term    = 5.61d+6 * t9a56 * tf%t9i32 * exp(-12.826/t9a13)
    dtermdt = term*(dt9a56/t9a56 - 1.5d0*tf%t9i &
         + 12.826/t9a13**2 * dt9a13)

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.11e+10 * tf%t932 * exp(-18.423*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 18.423*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_he3he4



  subroutine rate_c12pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee

    double precision, parameter :: q1 = 1.0d0/2.25d0

    !$gpu

    ! c12(p,g)13n
    aa   = 2.04e+07 * tf%t9i23 * exp(-13.69*tf%t9i13 - tf%t92*q1)
    daa  = aa*(-twoth*tf%t9i + oneth*13.69*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb   = 1.0d0 + 0.03*tf%t913 + 1.19*tf%t923 + 0.254*tf%t9 &
         + 2.06*tf%t943 + 1.12*tf%t953
    dbb  = oneth*0.03*tf%t9i23 + twoth*1.19*tf%t9i13 + 0.254 &
         + fourth*2.06*tf%t913 + fiveth*1.12*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 1.08e+05 * tf%t9i32 * exp(-4.925*tf%t9i)
    ddd  = dd*(-1.5d0*tf%t9i + 4.925*tf%t9i2)

    ee   = 2.15e+05 * tf%t9i32 * exp(-18.179*tf%t9i)
    dee  = ee*(-1.5d0*tf%t9i + 18.179*tf%t9i2)

    term    = cc + dd + ee
    dtermdt = dcc + ddd + dee

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term

    rev      = 8.84e+09 * tf%t932 * exp(-22.553*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 22.553*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    drrdd = 0.0d0

  end subroutine rate_c12pg



  subroutine rate_n14pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee

    double precision, parameter :: q1 = 1.0d0/10.850436d0

    !$gpu

    ! n14(p,g)o15
    aa  = 4.90e+07 * tf%t9i23 * exp(-15.228*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*15.228*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb   = 1.0d0 + 0.027*tf%t913 - 0.778*tf%t923 - 0.149*tf%t9 &
         + 0.261*tf%t943 + 0.127*tf%t953
    dbb  = oneth*0.027*tf%t9i23 - twoth*0.778*tf%t9i13 - 0.149 &
         + fourth*0.261*tf%t913 + fiveth*0.127*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 2.37e+03 * tf%t9i32 * exp(-3.011*tf%t9i)
    ddd  = dd*(-1.5d0*tf%t9i + 3.011*tf%t9i2)

    ee   = 2.19e+04 * exp(-12.530*tf%t9i)
    dee  = ee*12.530*tf%t9i2

    term    = cc + dd + ee
    dtermdt = dcc + ddd + dee

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev    = 2.70e+10 * tf%t932 * exp(-84.678*tf%t9i)
    drevdt = rev*(1.5d0*tf%t9i + 84.678*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_n14pg



  subroutine rate_n15pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,ff,dff

    double precision, parameter :: q1 = 1.0d0/0.2025d0

    !$gpu

    ! n15(p,g)o16
    aa  = 9.78e+08 * tf%t9i23 * exp(-15.251*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*15.251*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb   = 1.0d0  + 0.027*tf%t913 + 0.219*tf%t923 + 0.042*tf%t9 &
         + 6.83*tf%t943 + 3.32*tf%t953
    dbb  = oneth*0.027*tf%t9i23 + twoth*0.219*tf%t9i13 + 0.042 &
         + fourth*6.83*tf%t913 + fiveth*3.32*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 1.11e+04*tf%t9i32*exp(-3.328*tf%t9i)
    ddd  = dd*(-1.5d0*tf%t9i + 3.328*tf%t9i2)

    ee   = 1.49e+04*tf%t9i32*exp(-4.665*tf%t9i)
    dee  = ee*(-1.5d0*tf%t9i + 4.665*tf%t9i2)

    ff   = 3.8e+06*tf%t9i32*exp(-11.048*tf%t9i)
    dff  = ff*(-1.5d0*tf%t9i + 11.048*tf%t9i2)

    term    = cc + dd + ee + ff
    dtermdt = dcc + ddd + dee + dff

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 3.62e+10 * tf%t932 * exp(-140.734*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 140.734*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_n15pg



  subroutine rate_n15pa(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg

    double precision, parameter :: theta = 0.1d0
    double precision, parameter :: q1    = 1.0d0/0.272484d0

    !$gpu

    ! n15(p,a)c12
    aa  = 1.08d+12*tf%t9i23*exp(-15.251*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*15.251*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb   = 1.0d0 + 0.027*tf%t913 + 2.62*tf%t923 + 0.501*tf%t9 &
         + 5.36*tf%t943 + 2.60*tf%t953
    dbb  = oneth*0.027*tf%t9i23 + twoth*2.62*tf%t9i13 + 0.501 &
         + fourth*5.36*tf%t913 + fiveth*2.60*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 1.19d+08 * tf%t9i32 * exp(-3.676*tf%t9i)
    ddd  = dd*(-1.5d0*tf%t9i + 3.676*tf%t9i2)

    ee   = 5.41d+08 * tf%t9i12 * exp(-8.926*tf%t9i)
    dee  = ee*(-0.5d0*tf%t9i + 8.926*tf%t9i2)

    ff   = theta * 4.72d+08 * tf%t9i32 * exp(-7.721*tf%t9i)
    dff  = ff*(-1.5d0*tf%t9i + 7.721*tf%t9i2)

    gg   = theta * 2.20d+09 * tf%t9i32 * exp(-11.418*tf%t9i)
    dgg  = gg*(-1.5d0*tf%t9i + 11.418*tf%t9i2)

    term    = cc + dd + ee + ff + gg
    dtermdt = dcc + ddd + dee + dff + dgg

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 7.06d-01*exp(-57.625*tf%t9i)
    drevdt   = rev*57.625*tf%t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = rev * term

  end subroutine rate_n15pa



  subroutine rate_o16pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,zz

    !$gpu

    ! o16(p,g)f17
    aa  = exp(-0.728*tf%t923)
    daa = -twoth*aa*0.728*tf%t9i13

    bb  = 1.0d0 + 2.13 * (1.0d0 - aa)
    dbb = -2.13*daa

    cc  = tf%t923 * bb
    dcc = twoth*cc*tf%t9i + tf%t923*dbb

    dd   = exp(-16.692*tf%t9i13)
    ddd  = oneth*dd*16.692*tf%t9i43

    zz   = 1.0d0/cc
    ee   = dd*zz
    dee  = (ddd - ee*dcc)*zz

    term    = 1.50d+08 * ee
    dtermdt = 1.50d+08 * dee


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 3.03e+09*tf%t932*exp(-6.968*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 6.968*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_o16pg



  subroutine rate_n14ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,ff,dff

    double precision, parameter :: q1 = 1.0d0/0.776161d0

    !$gpu

    ! n14(a,g)f18
    aa  = 7.78d+09 * tf%t9i23 * exp(-36.031*tf%t9i13- tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*36.031*tf%t9i43 - 2.0d0*tf%t9*q1)

    bb   = 1.0d0 + 0.012*tf%t913 + 1.45*tf%t923 + 0.117*tf%t9 &
         + 1.97*tf%t943 + 0.406*tf%t953
    dbb  = oneth*0.012*tf%t9i23 + twoth*1.45*tf%t9i13 + 0.117 &
         + fourth*1.97*tf%t913 + fiveth*0.406*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 2.36d-10 * tf%t9i32 * exp(-2.798*tf%t9i)
    ddd  = dd*(-1.5d0*tf%t9i + 2.798*tf%t9i2)

    ee   = 2.03 * tf%t9i32 * exp(-5.054*tf%t9i)
    dee  = ee*(-1.5d0*tf%t9i + 5.054*tf%t9i2)

    ff   = 1.15d+04 * tf%t9i23 * exp(-12.310*tf%t9i)
    dff  = ff*(-twoth*tf%t9i + 12.310*tf%t9i2)

    term    = cc + dd + ee + ff
    dtermdt = dcc + ddd + dee + dff

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 5.42e+10 * tf%t932 * exp(-51.236*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 51.236*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_n14ag



  subroutine rate_fe52ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,tq2

    !$gpu

    ! fe52(n,g)fe53
    tq2     = tf%t9 - 0.348d0
    term    = 9.604d+05 * exp(-0.0626*tq2)
    dtermdt = -term*0.0626

    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 2.43d+09 * tf%t932 * exp(-123.951*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 123.951*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_fe52ng



  subroutine rate_fe53ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,tq1,tq10,dtq10,tq2

    !$gpu

    ! fe53(n,g)fe54
    tq1   = tf%t9/0.348
    tq10  = tq1**0.10
    dtq10 = 0.1d0*tq10/(0.348*tq1)
    tq2   = tf%t9 - 0.348d0

    term    = 1.817d+06 * tq10 * exp(-0.06319*tq2)
    dtermdt = term/tq10*dtq10 - term*0.06319

    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 1.56d+11 * tf%t932 * exp(-155.284*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 155.284*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_fe53ng



  subroutine rate_fe54ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den, fr, dfrdt, dfrdd, rr, drrdt, drrdd
    type (tf_t)      :: tf

    double precision :: aa, daa, bb, dbb, term, dtermdt

    !$gpu

    ! fe54(n,g)fe55
    aa   =  2.307390d+01 - 7.931795d-02 * tf%t9i + 7.535681d+00 * tf%t9i13 &
         - 1.595025d+01 * tf%t913 + 1.377715d+00 * tf%t9 - 1.291479d-01 * tf%t953 &
         + 6.707473d+00 * log(tf%t9)

    daa  =  7.931795d-02 * tf%t9i2 - oneth * 7.535681d+00 * tf%t9i43 &
         - oneth * 1.595025d+01 *tf%t9i23 + 1.377715d+00 - fiveth * 1.291479d-01 *tf%t923 &
         + 6.707473d+00 * tf%t9i

    if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
    else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
    end if

    bb  = 4.800293d+09 * tf%t932 * exp(-1.078986d+02 * tf%t9i)
    dbb = bb*(1.5d0*tf%t9i + 1.078986d+02 * tf%t9i2)

    ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0d-9
    !drrdd = 0.0d0

    ! forward rate
    !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den

  end subroutine rate_fe54ng



  subroutine rate_fe54pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: term,dtermdt,rev,drevdt,aa,daa,z,z2,z3

    !$gpu

    ! fe54(p,g)co55
    z     = min(tf%t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 9.593d-2*z - 3.445d-3*z2 + 8.594d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 9.593d-2 - 2.0d0*3.445d-3*tf%t9 + 3.0d0*8.594d-5*tf%t92
    end if

    term    = 4.51d+17 * tf%t9i23 * exp(-38.483 * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 38.483*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    !dfrdd = term

    rev      = 2.400d+09 * tf%t932 * exp(-58.605*tf%t9i)
    drevdt   = rev*(1.5d0*tf%t9i + 58.605*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    !drrdd = 0.0d0

  end subroutine rate_fe54pg




  subroutine rate_fe54ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: aa,daa,bb,dbb,term,dtermdt

    !$gpu

    ! fe54(a,p)co57
    aa   =  3.97474900d+01 - 6.06543100d+00 * tf%t9i + 1.63239600d+02 * tf%t9i13 &
         - 2.20457700d+02 * tf%t913 + 8.63980400d+00 * tf%t9 - 3.45841300d-01 * tf%t953 &
         + 1.31464200d+02 * log(tf%t9)

    daa  =  6.06543100d+00 * tf%t9i2 - oneth * 1.63239600d+02 * tf%t9i43 &
         - oneth * 2.20457700d+02 * tf%t9i23 + 8.63980400d+00 - fiveth * 3.45841300d-01 * tf%t923 &
         + 1.31464200d+02  * tf%t9i

    if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
    else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
    end if

    bb  = 2.16896000d+00  * exp(-2.05631700d+01 * tf%t9i)
    dbb = bb * 2.05631700d+01 * tf%t9i2

    ! reverse rate
    !drrdd = term
    rr    = term*den
    drrdt = dtermdt*den

    ! forward rate
    fr    = rr*bb
    dfrdt = drrdt*bb + rr*dbb*1.0d-9
    !dfrdd = drrdd*bb

  end subroutine rate_fe54ap



  subroutine rate_fe55ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: aa,daa,bb,dbb,term,dtermdt

    !$gpu

    ! fe55(n,g)fe56
    aa   =  1.954115d+01 - 6.834029d-02 * tf%t9i + 5.379859d+00 * tf%t9i13 &
         - 8.758150d+00 * tf%t913 + 5.285107d-01 * tf%t9 - 4.973739d-02  * tf%t953 &
         + 4.065564d+00  * log(tf%t9)

    daa  =  6.834029d-02 * tf%t9i2 - oneth * 5.379859d+00 * tf%t9i43 &
         - oneth * 8.758150d+00 * tf%t9i23 + 5.285107d-01 - fiveth * 4.973739d-02  *tf%t923 &
         + 4.065564d+00  * tf%t9i

    if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
    else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
    end if

    bb  = 7.684279d+10  * tf%t932 * exp(-1.299472d+02  * tf%t9i)
    dbb = bb*(1.5d0*tf%t9i + 1.299472d+02 * tf%t9i2)

    ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0d-9
    !drrdd = 0.0d0

    ! forward rate
    !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den

  end subroutine rate_fe55ng




  subroutine rate_fe56pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    double precision :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    double precision :: aa,daa,bb,dbb,term,dtermdt

    !$gpu

    ! fe56(p,g)co57

    aa   =  1.755960d+02 - 7.018872d+00 * tf%t9i + 2.800131d+02 * tf%t9i13 &
         - 4.749343d+02 * tf%t913 + 2.683860d+01 * tf%t9 - 1.542324d+00  * tf%t953 &
         + 2.315911d+02  * log(tf%t9)

    daa  =  7.018872d+00 * tf%t9i2 - oneth * 2.800131d+02 * tf%t9i43 &
         - oneth * 4.749343d+02 * tf%t9i23 + 2.683860d+01 - fiveth * 1.542324d+00  *tf%t923 &
         + 2.315911d+02  * tf%t9i

    if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
    else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
    end if

    bb  = 2.402486d+09 * tf%t932 * exp(-6.995192d+01 * tf%t9i)
    dbb = bb*(1.5d0*tf%t9i + 6.995192d+01 * tf%t9i2)


    ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0d-9
    !drrdd = 0.0d0

    ! forward rate
    !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den

  end subroutine rate_fe56pg



  ! this routine evaluates Langanke et al. 2000 fits for the ni56 electron
  ! capture rate rn56ec and neutrino loss rate sn56ec

  ! input:
  ! y56 = nickel56 molar abundance
  ! ye  = electron to baryon number, zbar/abar

  ! output:
  ! rn56ec = ni56 electron capture rate
  ! sn56ec = ni56 neutrino loss rate

  subroutine langanke(btemp,bden,y56,ye,rn56ec,sn56ec)

    implicit none

    integer          :: jp,kp,jr,jd
    double precision :: btemp,bden,y56,ye,rn56ec,sn56ec

    double precision :: rnt(2),rne(2,14),t9,r,rfm,rf0, &
                        rf1,rf2,dfacm,dfac0,dfac1,dfac2, &
                        tfm,tf0,tf1,tf2,tfacm,tfac0,tfac1,tfac2

    !$gpu

    ! calculate ni56 electron capture and neutrino loss rates
    rn56ec = 0.0
    sn56ec = 0.0
    if ( (btemp .lt. 1.0e9) .or. (bden*ye .lt. 1.0e6)) return
    t9    = min(btemp,1.4d10) * 1.0d-9
    r     = max(6.0d0,min(11.0d0,log10(bden*ye)))
    jp    = min(max(2,int(t9)),12)
    kp    = min(max(2,int(r)-5),4)
    rfm   = r - rv(kp-1)
    rf0   = r - rv(kp)
    rf1   = r - rv(kp+1)
    rf2   = r - rv(kp+2)
    dfacm = rf0*rf1*rf2*rfdm(kp)
    dfac0 = rfm*rf1*rf2*rfd0(kp)
    dfac1 = rfm*rf0*rf2*rfd1(kp)
    dfac2 = rfm*rf0*rf1*rfd2(kp)
    tfm   = t9 - tv(jp-1)
    tf0   = t9 - tv(jp)
    tf1   = t9 - tv(jp+1)
    tf2   = t9 - tv(jp+2)
    tfacm = tf0*tf1*tf2*tfdm(jp)
    tfac0 = tfm*tf1*tf2*tfd0(jp)
    tfac1 = tfm*tf0*tf2*tfd1(jp)
    tfac2 = tfm*tf0*tf1*tfd2(jp)

    ! evaluate the spline fits
    do jr = 1,2
       do jd = jp-1,jp+2
          rne(jr,jd) =   dfacm*datn(jr,kp-1,jd) + dfac0*datn(jr,kp,jd) &
               + dfac1*datn(jr,kp+1,jd) + dfac2*datn(jr,kp+2,jd)
       enddo
       rnt(jr) =  tfacm*rne(jr,jp-1) + tfac0*rne(jr,jp) &
            + tfac1*rne(jr,jp+1) + tfac2*rne(jr,jp+2)
    enddo

    ! set the output
    rn56ec = 10.0d0**rnt(1)
    sn56ec = 6.022548d+23 * 1.60218d-6 * y56 * 10.0d0**rnt(2)
   
    !write(*,*) "btemp",btemp, "bden", bden, "t9",t9,"r",r, "rn56ec",rn56ec, "sn56ec", sn56ec

  end subroutine langanke



  ! given the electron degeneracy parameter etakep (chemical potential
  ! without the electron's rest mass divided by kt) and the temperature
  ! temp, this routine calculates rates for
  ! electron capture on protons rpen (captures/sec/proton),
  ! positron capture on neutrons rnep (captures/sec/neutron),
  ! and their associated neutrino energy loss rates
  ! spenc (erg/sec/proton) and snepc (erg/sec/neutron)

  subroutine ecapnuc(etakep,temp,rpen,rnep,spenc,snepc)

    implicit none

    double precision :: etakep,temp,rpen,rnep,spenc,snepc

    integer          iflag
    double precision :: t9,t5,qn,etaef,etael,zetan,eta,etael2, &
                        etael3,etael4,f1l,f2l,f3l,f4l,f5l,f1g, &
                        f2g,f3g,f4g,f5g,exmeta,eta2,eta3,eta4, &
                        fac0,fac1,fac2,fac3,rie1,rie2,facv0,facv1, &
                        facv2,facv3,facv4,rjv1,rjv2,spen,snep, &
                        pi2,exeta,zetan2,f0,etael5,bktinv, &
                        qn1,ftinv,twoln,cmk5,cmk6,bk,pi,qn2,c2me, &
                        xmp,xmn,qndeca,tmean
    parameter        (qn1    = -2.0716446d-06, &
         ftinv  = 1.0d0/1083.9269d0, &
         twoln  = 0.6931472d0, &
         cmk5   = 1.3635675d-49, &
         cmk6   = 2.2993864d-59, &
         bk     = 1.38062e-16, &
         pi     = 3.1415927d0, &
         pi2    = pi * pi, &
         qn2    = 2.0716446d-06, &
         c2me   = 8.1872665d-07, &
         xmp    = 1.6726485d-24, &
         xmn    = 1.6749543d-24, &
         qndeca = 1.2533036d-06, &
         tmean  = 886.7d0)
    !     3                  tmean  = 935.14d0)

    double precision third,sixth
    parameter        (third = 1.0d0/3.0d0, &
         sixth = 1.0d0/6.0d0)

    !$gpu

    ! tmean and qndeca are the mean lifetime and decay energy of the neutron
    ! xmp,xnp are masses of the p and n in grams.
    ! c2me is the constant used to convert the neutrino energy
    ! loss rate from mec2/s (as in the paper) to ergs/particle/sec.

    ! initialize
    rpen   = 0.0d0
    rnep   = 0.0d0
    spen   = 0.0d0
    snep   = 0.0d0
    t9     = temp * 1.0d-9
    bktinv = 1.0d0/(bk *temp)
    iflag  = 0
    qn     = qn1


    ! chemical potential including the electron rest mass
    etaef = etakep + c2me*bktinv


    ! iflag=1 is for electrons,  iflag=2 is for positrons
502 iflag = iflag + 1
    if (iflag.eq.1) etael = qn2*bktinv
    if (iflag.eq.2) then
       etael = c2me*bktinv
       etaef = -etaef
    endif

    t5    = temp*temp*temp*temp*temp
    zetan = qn*bktinv
    eta   = etaef - etael

    ! protect from overflowing with large eta values
    if (eta .le. 6.8e+02) then
       exeta = exp(eta)
    else
       exeta = 0.0d0
    end if
    etael2 = etael*etael
    etael3 = etael2*etael
    etael4 = etael3*etael
    etael5 = etael4*etael
    zetan2 = zetan*zetan
    if (eta .le. 6.8e+02) then
       f0 = log(1.0d0 + exeta)
    else
       f0 = eta
    end if

    ! if eta le. 0., the following fermi integrals apply
    f1l = exeta
    f2l = 2.0d0   * f1l
    f3l = 6.0d0   * f1l
    f4l = 24.0d0  * f1l
    f5l = 120.0d0 * f1l

    ! if eta gt. 0., the following fermi integrals apply:
    f1g = 0.0d0
    f2g = 0.0d0
    f3g = 0.0d0
    f4g = 0.0d0
    f5g = 0.0d0
    if (eta .gt. 0.0) then
       exmeta = dexp(-eta)
       eta2   = eta*eta
       eta3   = eta2*eta
       eta4   = eta3*eta
       f1g = 0.5d0*eta2 + 2.0d0 - exmeta
       f2g = eta3*third + 4.0d0*eta + 2.0d0*exmeta
       f3g = 0.25d0*eta4 + 0.5d0*pi2*eta2 + 12.0d0 - 6.0d0*exmeta
       f4g = 0.2d0*eta4*eta + 2.0d0*pi2*third*eta3 + 48.0d0*eta &
            + 24.0d0*exmeta
       f5g = eta4*eta2*sixth + 5.0d0*sixth*pi2*eta4 &
            + 7.0d0*sixth*pi2*eta2  + 240.0d0 -120.d0*exmeta
    end if

    ! factors which are multiplied by the fermi integrals
    fac3 = 2.0d0*zetan + 4.0d0*etael
    fac2 = 6.0d0*etael2 + 6.0d0*etael*zetan + zetan2
    fac1 = 4.0d0*etael3 + 6.0d0*etael2*zetan + 2.0d0*etael*zetan2
    fac0 = etael4 + 2.0d0*zetan*etael3 + etael2*zetan2

    ! electron capture rates onto protons with no blocking
    rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
    rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

    ! neutrino emission rate for electron capture:
    facv4 = 5.0d0*etael + 3.0d0*zetan
    facv3 = 10.0d0*etael2 + 12.0d0*etael*zetan + 3.0d0*zetan2
    facv2 = 10.0d0*etael3 + 18.0d0*etael2*zetan &
         + 9.0d0*etael*zetan2 + zetan2*zetan
    facv1 = 5.0d0*etael4 + 12.0d0*etael3*zetan &
         + 9.0d0*etael2*zetan2 + 2.0d0*etael*zetan2*zetan
    facv0 = etael5 + 3.0d0*etael4*zetan &
         + 3.0d0*etael3*zetan2 + etael2*zetan2*zetan
    rjv1  = f5l + facv4*f4l + facv3*f3l &
         + facv2*f2l + facv1*f1l + facv0*f0
    rjv2  = f5g + facv4*f4g + facv3*f3g &
         + facv2*f2g + facv1*f1g + facv0*f0

    ! for electrons capture onto protons
    if (iflag.eq.2) go to 503
    if (eta.gt.0.) go to 505
    rpen  = twoln*cmk5*t5*rie1*ftinv
    spen  = twoln*cmk6*t5*temp*rjv1*ftinv
    spenc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me
    go to 504
505 rpen = twoln*cmk5*t5*rie2*ftinv
    spen = twoln*cmk6*t5*temp*rjv2*ftinv
    spenc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me
504 continue
    qn = qn2
    go to 502

    ! for positrons capture onto neutrons
503 if (eta.gt.0.) go to 507
    rnep  = twoln*cmk5*t5*rie1*ftinv
    snep  = twoln*cmk6*t5*temp*rjv1*ftinv
    snepc = twoln*cmk6*t5*temp*rjv1*ftinv*c2me
    go to 506
507 rnep  = twoln*cmk5*t5*rie2*ftinv
    snep  = twoln*cmk6*t5*temp*rjv2*ftinv
    snepc = twoln*cmk6*t5*temp*rjv2*ftinv*c2me
506 continue
    return
  end subroutine ecapnuc

end module aprox_rates_module
