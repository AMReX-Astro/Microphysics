module aprox_rates_module

  ! These rate routines come from the public_aprox13/19/21.f90 files.
  ! Only those used in the networks we have in this repository are kept.
  ! We modify the calling sequence to take the tfactors as an argument.
  ! We comment out the density derivatives because we never evolve density
  ! in our reaction networks.

  use tfactors_module
  use microphysics_type_module, only: rt

  implicit none

  real(rt), allocatable :: rv(:), tv(:), datn(:,:,:)
  real(rt), allocatable :: rfdm(:),rfd0(:),rfd1(:),rfd2(:)
  real(rt), allocatable :: tfdm(:),tfd0(:),tfd1(:),tfd2(:)

#ifdef AMREX_USE_CUDA
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

    rv = (/ 6.0_rt, 7.0_rt, 8.0_rt, 9.0_rt, 10.0_rt, 11.0_rt /)
    tv = (/ 1.0_rt,2.0_rt,3.0_rt,4.0_rt,5.0_rt,6.0_rt,7.0_rt,8.0_rt,9.0_rt,10.0_rt,11.0_rt,12.0_rt,13.0_rt,14.0_rt /)

    datn(1,:,:) = reshape( (/ -4.363_rt, -3.091_rt, -1.275_rt, 1.073_rt, 3.035_rt, 4.825_rt, &
			      -4.17_rt, -2.964_rt, -1.177_rt, 1.085_rt, 3.037_rt, 4.826_rt, &
			      -3.834_rt, -2.727_rt, -1.039_rt, 1.104_rt, 3.04_rt, 4.826_rt, &
			      -3.284_rt, -2.418_rt, -0.882_rt, 1.129_rt, 3.043_rt, 4.827_rt, &
			      -2.691_rt, -2.093_rt, -0.719_rt, 1.159_rt, 3.048_rt, 4.827_rt, &
			      -2.1675_rt, -1.7668_rt, -0.5573_rt, 1.1947_rt, 3.0527_rt, 4.8272_rt, &
			      -1.7095_rt, -1.4462_rt, -0.3991_rt, 1.2358_rt, 3.0577_rt, 4.8276_rt, &
        		      -1.3119_rt, -1.1451_rt, -0.2495_rt, 1.2818_rt, 3.0648_rt, 4.8284_rt, &
        		      -0.9812_rt, -0.8612_rt, -0.1084_rt, 1.3336_rt, 3.0738_rt, 4.8295_rt, &
        		      -0.682_rt, -0.595_rt, 0.028_rt, 1.386_rt, 3.084_rt, 4.831_rt, &
                              -0.4046_rt, -0.3523_rt, 0.1605_rt, 1.4364_rt, 3.0957_rt, 4.8333_rt, &
                              -0.1636_rt, -0.1352_rt, 0.2879_rt, 1.4861_rt, 3.1092_rt, 4.8365_rt, &
                              0.0461_rt, 0.0595_rt, 0.4105_rt, 1.5354_rt, 3.1242_rt, 4.8405_rt, &
                              0.2295_rt, 0.235_rt, 0.5289_rt, 1.5842_rt, 3.1405_rt, 4.845_rt /), &
                           (/ 6, 14 /) )

    datn(2,:,:) = reshape( (/ -4.539_rt, -3.097_rt, -1.134_rt, 1.525_rt, 3.907_rt, 6.078_rt, &
        		      -4.199_rt, -2.905_rt, -1.024_rt, 1.545_rt, 3.91_rt, 6.079_rt, &
        		      -3.736_rt, -2.602_rt, -0.851_rt, 1.578_rt, 3.916_rt, 6.08_rt, &
        		      -3.052_rt, -2.206_rt, -0.636_rt, 1.623_rt, 3.923_rt, 6.081_rt, &
        		      -2.31_rt, -1.766_rt, -0.396_rt, 1.678_rt, 3.931_rt, 6.082_rt, &
        		      -1.6631_rt, -1.319_rt, -0.1438_rt, 1.7471_rt, 3.9409_rt, 6.0829_rt, &
        		      -1.1064_rt, -0.8828_rt, 0.1094_rt, 1.8279_rt, 3.9534_rt, 6.0841_rt, &
       			      -0.6344_rt, -0.496_rt, 0.3395_rt, 1.9168_rt, 3.9699_rt, 6.0862_rt, &
        		      -0.2568_rt, -0.1555_rt, 0.5489_rt, 2.0163_rt, 3.9906_rt, 6.0893_rt, &
        		      0.081_rt, 0.158_rt, 0.746_rt, 2.114_rt, 4.013_rt, 6.093_rt, &
        		      0.3961_rt, 0.4448_rt, 0.9304_rt, 2.2026_rt, 4.0363_rt, 6.0976_rt, &
        		      0.6673_rt, 0.6964_rt, 1.0985_rt, 2.2849_rt, 4.0614_rt, 6.1033_rt, &
        		      0.9009_rt, 0.9175_rt, 1.2525_rt, 2.3619_rt, 4.0882_rt, 6.1099_rt, &
        		      1.1032_rt, 1.113_rt, 1.3947_rt, 2.4345_rt, 4.1161_rt, 6.1171_rt /), &
                           (/ 6, 14 /) )



    ! Evaluate the cubic interp parameters for ni56 electron capture
    ! which is used in the langanke subroutine.

    do k = 2, 4
       rfdm(k)=1._rt/((rv(k-1)-rv(k))*(rv(k-1)-rv(k+1))*(rv(k-1)-rv(k+2)))
       rfd0(k)=1._rt/((rv(k)-rv(k-1))*(rv(k)-rv(k+1))*(rv(k)-rv(k+2)))
       rfd1(k)=1._rt/((rv(k+1)-rv(k-1))*(rv(k+1)-rv(k))*(rv(k+1)-rv(k+2)))
       rfd2(k)=1._rt/((rv(k+2)-rv(k-1))*(rv(k+2)-rv(k))*(rv(k+2)-rv(k+1)))
    enddo

    do j = 2, 12
       tfdm(j)=1._rt/((tv(j-1)-tv(j))*(tv(j-1)-tv(j+1))*(tv(j-1)-tv(j+2)))
       tfd0(j)=1._rt/((tv(j)-tv(j-1))*(tv(j)-tv(j+1))*(tv(j)-tv(j+2)))
       tfd1(j)=1._rt/((tv(j+1)-tv(j-1))*(tv(j+1)-tv(j))*(tv(j+1)-tv(j+2)))
       tfd2(j)=1._rt/((tv(j+2)-tv(j-1))*(tv(j+2)-tv(j))*(tv(j+2)-tv(j+1)))
    enddo

  !$acc update device(rv, tv, datn, rfdm, rfd0, rfd1, rfd2, tfdm, tfd0, tfd1, tfd2)

  end subroutine rates_init



  subroutine rate_c12ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                        dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2, &
                        zz

    real(rt), parameter :: q1 = 1.0e0_rt/12.222016e0_rt

    !$gpu

    ! c12(a,g)o16
    aa   = 1.0e0_rt + 0.0489e0_rt*tf%t9i23
    daa  = -twoth*0.0489e0_rt*tf%t9i53

    bb   = tf%t92*aa*aa
    dbb  = 2.0e0_rt*(bb*tf%t9i + tf%t92*aa*daa)

    cc   = exp(-32.120e0_rt*tf%t9i13 - tf%t92*q1)
    dcc  = cc * (oneth*32.120e0_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    dd   = 1.0e0_rt + 0.2654e0_rt*tf%t9i23
    ddd  = -twoth*0.2654e0_rt*tf%t9i53

    ee   = tf%t92*dd*dd
    dee  = 2.0e0_rt*(ee*tf%t9i + tf%t92*dd*ddd)

    ff   = exp(-32.120e0_rt*tf%t9i13)
    dff  = ff * oneth*32.120e0_rt*tf%t9i43

    gg   = 1.25e3_rt * tf%t9i32 * exp(-27.499_rt*tf%t9i)
    dgg  = gg*(-1.5e0_rt*tf%t9i + 27.499_rt*tf%t9i2)

    hh   = 1.43e-2_rt * tf%t95 * exp(-15.541_rt*tf%t9i)
    dhh  = hh*(5.0e0_rt*tf%t9i + 15.541_rt*tf%t9i2)

    zz   = 1.0e0_rt/bb
    f1   = cc*zz
    df1  = (dcc - f1*dbb)*zz

    zz   = 1.0e0_rt/ee
    f2   = ff*zz
    df2  = (dff - f2*dee)*zz

    term    = 1.04e8_rt*f1  + 1.76e8_rt*f2 + gg + hh
    dtermdt = 1.04e8_rt*df1 + 1.76e8_rt*df2 + dgg + dhh


    ! 1.7_rt times cf88 value
    term     = 1.7e0_rt * term
    dtermdt  = 1.7e0_rt * dtermdt

    fr    = term * den
    dfrdt = dtermdt * den * 1.0e-9_rt
    !dfrdd = term

    rev    = 5.13e10_rt * tf%t932 * exp(-83.111_rt*tf%t9i)
    drevdt = rev*(1.5e0_rt*tf%t9i + 83.111_rt*tf%t9i2)

    rr     = rev * term
    drrdt  = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd  = 0.0e0_rt

  end subroutine rate_c12ag

  ! This routine computes the nuclear reaction rate for 12C(a,g)16O and its inverse 
  ! using fit parameters from Deboer et al. 2017 (https://doi.org/10.1103_rt/RevModPhys.89.035007_rt).
  subroutine rate_c12ag_deboer17(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    type (tf_t)      :: tf

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd, &
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
    a0_nr = 24.1e0_rt
    a1_nr = 0.e0_rt 
    a2_nr = -32.e0_rt
    a3_nr = -5.9e0_rt
    a4_nr = 1.8e0_rt
    a5_nr = -0.17e0_rt
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

    dterm_a0_nr = 0.e0_rt
    dterm_a1_nr = 0.e0_rt
    dterm_a2_nr = -a2_nr*tf%t9i43*term_a2_nr/3.e0_rt
    dterm_a3_nr = a3_nr*tf%t9i23*term_a3_nr/3.e0_rt
    dterm_a4_nr = a4_nr*term_a4_nr
    dterm_a5_nr = a5_nr*tf%t923*term_a5_nr*fiveth
    dterm_a6_nr = tf%t9i*a6_nr*tf%t9**a6_nr

    dterm_nr = (term_a0_nr * term_a1_nr * dterm_a2_nr * term_a3_nr * term_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * dterm_a3_nr * term_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * dterm_a4_nr * term_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * term_a4_nr * dterm_a5_nr * term_a6_nr) + &
               (term_a0_nr * term_a1_nr * term_a2_nr * term_a3_nr * term_a4_nr * term_a5_nr * dterm_a6_nr)

    ! resonant contributions to the reaction
    a0_r = 7.4e0_rt
    a1_r = -30.e0_rt
    a2_r = 0.e0_rt
    a3_r = 0.e0_rt
    a4_r = 0.e0_rt
    a5_r = 0.e0_rt
    a6_r = -3.0e0_rt/2.0e0_rt

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

    dterm_a0_r = 0.e0_rt
    dterm_a1_r = -a1_r*tf%t9i2*term_a1_r
    dterm_a2_r = 0.e0_rt
    dterm_a3_r = 0.e0_rt
    dterm_a4_r = 0.e0_rt
    dterm_a5_r = 0.e0_rt
    dterm_a6_r = tf%t9i*a6_r*tf%t9**a6_r
    
    dterm_r = (term_a0_r * dterm_a1_r * term_a6_r) + &
	      (term_a0_r * term_a1_r * dterm_a6_r)
    

    ! full rate is the sum of resonant and non-resonant contributions
    term = term_nr + term_r 
    dtermdt = dterm_nr + dterm_r

    fr    = term * den
    dfrdt = dtermdt * den * 1.0e-9_rt

    ! first term is 9.8685e9_rt * T9**(2/3) * (M0*M1/M3)**(3/2) 
    ! see iliadis 2007 eqn. 3.44
    ! ratio of partition functions are assumed to be unity
    rev    = 5.1345573e10_rt * tf%t932 * exp(-83.114082_rt*tf%t9i)
    drevdt = rev*(1.5e0_rt*tf%t9i + 83.114082_rt*tf%t9i2)

    rr     = rev * term
    drrdt  = (drevdt*term + rev*dtermdt) * 1.0e-9_rt


  end subroutine rate_c12ag_deboer17


  subroutine rate_tripalf(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac, &
                        drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1

    real(rt), parameter :: rc28   = 0.1e0_rt
    real(rt), parameter :: q1     = 1.0e0_rt/0.009604e0_rt
    real(rt), parameter :: q2     = 1.0e0_rt/0.055225e0_rt

    !$gpu

    ! triple alfa to c12
    ! this is a(a,g)be8
    aa    = 7.40e05_rt * tf%t9i32 * exp(-1.0663_rt*tf%t9i)
    daa   = aa*(-1.5e0_rt*tf%t9i  + 1.0663_rt*tf%t9i2)

    bb    = 4.164e09_rt * tf%t9i23 * exp(-13.49_rt*tf%t9i13 - tf%t92*q1)
    dbb   = bb*(-twoth*tf%t9i + oneth*13.49_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    cc    = 1.0e0_rt + 0.031_rt*tf%t913 + 8.009_rt*tf%t923 + 1.732_rt*tf%t9 &
          + 49.883_rt*tf%t943 + 27.426_rt*tf%t953
    dcc   = oneth*0.031_rt*tf%t9i23 + twoth*8.009_rt*tf%t9i13 + 1.732_rt &
          + fourth*49.883_rt*tf%t913 + fiveth*27.426_rt*tf%t923

    r2abe    = aa + bb * cc
    dr2abedt = daa + dbb*cc + bb*dcc


    ! this is be8(a,g)c12
    dd    = 130.0e0_rt * tf%t9i32 * exp(-3.3364_rt*tf%t9i)
    ddd   = dd*(-1.5e0_rt*tf%t9i + 3.3364_rt*tf%t9i2)

    ee    = 2.510e07_rt * tf%t9i23 * exp(-23.57_rt*tf%t9i13 - tf%t92*q2)
    dee   = ee*(-twoth*tf%t9i + oneth*23.57_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q2)

    ff    = 1.0e0_rt + 0.018_rt*tf%t913 + 5.249_rt*tf%t923 + 0.650_rt*tf%t9 + &
         19.176_rt*tf%t943 + 6.034_rt*tf%t953
    dff   = oneth*0.018_rt*tf%t9i23 + twoth*5.249_rt*tf%t9i13 + 0.650_rt &
          + fourth*19.176_rt*tf%t913 + fiveth*6.034_rt*tf%t923

    rbeac    = dd + ee * ff
    drbeacdt = ddd + dee * ff + ee * dff


    ! a factor
    xx    = rc28 * 1.35e-07_rt * tf%t9i32 * exp(-24.811_rt*tf%t9i)
    dxx   = xx*(-1.5e0_rt*tf%t9i + 24.811_rt*tf%t9i2)


    ! high temperature rate
    if (tf%t9.gt.0.08_rt) then
       term    = 2.90e-16_rt * r2abe * rbeac + xx
       dtermdt =   2.90e-16_rt * dr2abedt * rbeac &
                 + 2.90e-16_rt * r2abe * drbeacdt &
                 + dxx


    ! low temperature rate
    else
       uu   = 0.8e0_rt*exp(-(0.025_rt*tf%t9i)**3.263_rt)
       yy   = 0.2e0_rt + uu
       ! fxt yy   = 0.01_rt + 0.2e0_rt + uu
       dyy  = uu * 3.263_rt*(0.025_rt*tf%t9i)**2.263_rt * (0.025_rt*tf%t9i2)
       vv   = 4.0e0_rt*exp(-(tf%t9/0.025_rt)**9.227_rt)
       zz   = 1.0e0_rt + vv
       dzz  = vv * 9.227_rt*(tf%t9/0.025_rt)**8.227_rt * 40.0e0_rt
       aa   = 1.0e0_rt/zz
       f1   = 0.01e0_rt + yy * aa
       ! fxt f1   = yy * aa
       df1  = (dyy - f1*dzz)*aa
       term = 2.90e-16_rt * r2abe * rbeac * f1 +  xx
       dtermdt =   2.90e-16_rt * dr2abedt * rbeac * f1 &
                 + 2.90e-16_rt * r2abe * drbeacdt * f1 &
                 + 2.90e-16_rt * r2abe * rbeac * df1 &
                 + dxx
    end if


    ! rates
    !      term    = 1.2e0_rt * term
    !      dtermdt = 1.2e0_rt * term

    fr    = term * den * den
    dfrdt = dtermdt * den * den * 1.0e-9_rt
    !dfrdd = 2.0e0_rt * term * den

    rev    = 2.00e20_rt*tf%t93*exp(-84.424_rt*tf%t9i)
    drevdt = rev*(3.0e0_rt*tf%t9i + 84.424_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_tripalf



  subroutine rate_c12c12(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                        aa,zz

    !$gpu

    ! c12 + c12 reaction
    aa      = 1.0e0_rt + 0.0396_rt*tf%t9
    zz      = 1.0e0_rt/aa

    t9a     = tf%t9*zz
    dt9a    = (1.0e0_rt -  t9a*0.0396_rt)*zz

    zz      = dt9a/t9a
    t9a13   = t9a**oneth
    dt9a13  = oneth*t9a13*zz

    t9a56   = t9a**fivsix
    dt9a56  = fivsix*t9a56*zz

    term    = 4.27e26_rt * t9a56 * tf%t9i32 * &
         exp(-84.165_rt/t9a13 - 2.12e-03_rt*tf%t93)
    dtermdt = term*(dt9a56/t9a56 - 1.5e0_rt*tf%t9i &
            + 84.165_rt/t9a13**2*dt9a13 - 6.36e-3_rt*tf%t92)

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_c12c12



  subroutine rate_c12o16(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23, &
                        t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz

    !$gpu

    ! c12 + o16 reaction; see cf88 references 47-4
    if (tf%t9.ge.0.5_rt) then
       aa     = 1.0e0_rt + 0.055_rt*tf%t9
       zz     = 1.0e0_rt/aa

       t9a    = tf%t9*zz
       dt9a   = (1.0e0_rt - t9a*0.055_rt)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0e0_rt * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa      = exp(-0.18_rt*t9a*t9a)
       daa     = -aa * 0.36_rt * t9a * dt9a

       bb      = 1.06e-03_rt*exp(2.562_rt*t9a23)
       dbb     = bb * 2.562_rt * dt9a23

       cc      = aa + bb
       dcc     = daa + dbb

       zz      = 1.0e0_rt/cc
       term    = 1.72e31_rt * t9a56 * tf%t9i32 * exp(-106.594_rt/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5e0_rt*tf%t9i &
                       + 106.594_rt/t9a23*dt9a13 - zz*dcc)

    else
       ! term    = 2.6288035e-29_rt
       term    = 0.0e0_rt
       dtermdt = 0.0e0_rt
    endif


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_c12o16



  subroutine rate_o16o16(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt

    !$gpu

    ! o16 + o16
    term  = 7.10e36_rt * tf%t9i23 * &
         exp(-135.93_rt * tf%t9i13 - 0.629_rt*tf%t923 &
         - 0.445_rt*tf%t943 + 0.0103_rt*tf%t9*tf%t9)

    dtermdt = -twoth*term*tf%t9i &
         + term * (oneth*135.93_rt*tf%t9i43 - twoth*0.629_rt*tf%t9i13 &
         - fourth*0.445_rt*tf%t913 + 0.0206_rt*tf%t9)


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_o16o16



  subroutine rate_o16ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                        cc,dcc,term2,dterm2,rev,drevdt

    real(rt), parameter :: q1 = 1.0e0_rt/2.515396e0_rt

    !$gpu

    ! o16(a,g)ne20
    term1   = 9.37e9_rt * tf%t9i23 * exp(-39.757_rt*tf%t9i13 - tf%t92*q1)
    dterm1  = term1*(-twoth*tf%t9i + oneth*39.757_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    aa      = 62.1_rt * tf%t9i32 * exp(-10.297_rt*tf%t9i)
    daa     = aa*(-1.5e0_rt*tf%t9i + 10.297_rt*tf%t9i2)

    bb      = 538.0e0_rt * tf%t9i32 * exp(-12.226_rt*tf%t9i)
    dbb     = bb*(-1.5e0_rt*tf%t9i + 12.226_rt*tf%t9i2)

    cc      = 13.0e0_rt * tf%t92 * exp(-20.093_rt*tf%t9i)
    dcc     = cc*(2.0e0_rt*tf%t9i + 20.093_rt*tf%t9i2)

    term2   = aa + bb + cc
    dterm2  = daa + dbb + dcc

    term    = term1 + term2
    dtermdt = dterm1 + dterm2


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 5.65e10_rt*tf%t932*exp(-54.937_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 54.937_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_o16ag



  subroutine rate_ne20ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                        term2,dterm2,term3,dterm3,rev,drevdt,zz

    real(rt), parameter :: rc102 = 0.1e0_rt
    real(rt), parameter :: q1    = 1.0e0_rt/4.923961e0_rt

    !$gpu

    ! ne20(a,g)mg24
    aa   = 4.11e11_rt * tf%t9i23 * exp(-46.766_rt*tf%t9i13 - tf%t92*q1)
    daa  = aa*(-twoth*tf%t9i + oneth*46.766_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb   = 1.0e0_rt + 0.009_rt*tf%t913 + 0.882_rt*tf%t923 + 0.055_rt*tf%t9 &
         + 0.749_rt*tf%t943 + 0.119_rt*tf%t953
    dbb  = oneth*0.009_rt*tf%t9i23 + twoth*0.882_rt*tf%t9i13 + 0.055_rt &
         + fourth*0.749_rt*tf%t913 + fiveth*0.119_rt*tf%t923

    term1  = aa * bb
    dterm1 = daa * bb + aa * dbb


    aa   = 5.27e03_rt * tf%t9i32 * exp(-15.869_rt*tf%t9i)
    daa  = aa*(-1.5e0_rt*tf%t9i + 15.869_rt*tf%t9i2)

    bb   = 6.51e03_rt * tf%t912 * exp(-16.223_rt*tf%t9i)
    dbb  = bb*(0.5e0_rt*tf%t9i + 16.223_rt*tf%t9i2)

    term2  = aa + bb
    dterm2 = daa + dbb


    aa   = 42.1_rt * tf%t9i32 * exp(-9.115_rt*tf%t9i)
    daa  = aa*(-1.5e0_rt*tf%t9i + 9.115_rt*tf%t9i2)

    bb   =  32.0_rt * tf%t9i23 * exp(-9.383_rt*tf%t9i)
    dbb  = bb*(-twoth*tf%t9i + 9.383_rt*tf%t9i2)

    term3  = rc102 * (aa + bb)
    dterm3 = rc102 * (daa + dbb)


    aa  = 5.0e0_rt*exp(-18.960_rt*tf%t9i)
    daa = aa*18.960_rt*tf%t9i2

    bb  = 1.0e0_rt + aa
    dbb = daa

    zz      = 1.0e0_rt/bb
    term    = (term1 + term2 + term3)*zz
    dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 6.01e10_rt * tf%t932 * exp(-108.059_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 108.059_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_ne20ag



  subroutine rate_mg24ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,gg,dgg,hh,hhi,rev,drevdt

    real(rt), parameter :: rc121 = 0.1e0_rt

    !$gpu

    ! 24mg(a,g)28si
    aa    = 4.78e01_rt * tf%t9i32 * exp(-13.506_rt*tf%t9i)
    daa   = aa*(-1.5e0_rt*tf%t9i + 13.506_rt*tf%t9i2)

    bb    =  2.38e03_rt * tf%t9i32 * exp(-15.218_rt*tf%t9i)
    dbb   = bb*(-1.5e0_rt*tf%t9i + 15.218_rt*tf%t9i2)

    cc    = 2.47e02_rt * tf%t932 * exp(-15.147_rt*tf%t9i)
    dcc   = cc*(1.5e0_rt*tf%t9i + 15.147_rt*tf%t9i2)

    dd    = rc121 * 1.72e-09_rt * tf%t9i32 * exp(-5.028_rt*tf%t9i)
    ddd   = dd*(-1.5e0_rt*tf%t9i + 5.028_rt*tf%t9i2)

    ee    = rc121* 1.25e-03_rt * tf%t9i32 * exp(-7.929_rt*tf%t9i)
    dee   = ee*(-1.5e0_rt*tf%t9i + 7.929_rt*tf%t9i2)

    ff    = rc121 * 2.43e01_rt * tf%t9i * exp(-11.523_rt*tf%t9i)
    dff   = ff*(-tf%t9i + 11.523_rt*tf%t9i2)

    gg    = 5.0e0_rt*exp(-15.882_rt*tf%t9i)
    dgg   = gg*15.882_rt*tf%t9i2

    hh    = 1.0e0_rt + gg
    hhi   = 1.0e0_rt/hh

    term    = (aa + bb + cc + dd + ee + ff) * hhi
    dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 6.27e10_rt * tf%t932 * exp(-115.862_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 115.862_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_mg24ag



  subroutine rate_mg24ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,gg,dgg,term1,dterm1,term2,dterm2, &
                        rev,drevdt

    real(rt), parameter :: rc148 = 0.1e0_rt
    real(rt), parameter :: q1    = 1.0e0_rt/0.024649e0_rt

    !$gpu

    ! 24mg(a,p)al27
    aa     = 1.10e08_rt * tf%t9i23 * exp(-23.261_rt*tf%t9i13 - tf%t92*q1)
    daa    = -twoth*aa*tf%t9i + aa*(23.261_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb     =  1.0e0_rt + 0.018_rt*tf%t913 + 12.85_rt*tf%t923 + 1.61_rt*tf%t9 &
         + 89.87_rt*tf%t943 + 28.66_rt*tf%t953
    dbb    = oneth*0.018_rt*tf%t9i23 + twoth*12.85_rt*tf%t9i13 + 1.61_rt &
           + fourth*89.87_rt*tf%t913 + fiveth*28.66_rt*tf%t923

    term1  = aa * bb
    dterm1 = daa * bb + aa * dbb

    aa     = 129.0e0_rt * tf%t9i32 * exp(-2.517_rt*tf%t9i)
    daa    = -1.5e0_rt*aa*tf%t9i + aa*2.517_rt*tf%t9i2

    bb     = 5660.0e0_rt * tf%t972 * exp(-3.421_rt*tf%t9i)
    dbb    = 3.5e0_rt*bb*tf%t9i +  bb*3.421_rt*tf%t9i2

    cc     = rc148 * 3.89e-08_rt * tf%t9i32 * exp(-0.853_rt*tf%t9i)
    dcc    = -1.5e0_rt*cc*tf%t9i + cc*0.853_rt*tf%t9i2

    dd     = rc148 * 8.18e-09_rt * tf%t9i32 * exp(-1.001_rt*tf%t9i)
    ddd    = -1.5e0_rt*dd*tf%t9i + dd*1.001_rt*tf%t9i2

    term2  = aa + bb + cc + dd
    dterm2 = daa + dbb + dcc + ddd

    ee     = oneth*exp(-9.792_rt*tf%t9i)
    dee    = ee*9.792_rt*tf%t9i2

    ff     =  twoth * exp(-11.773_rt*tf%t9i)
    dff    = ff*11.773_rt*tf%t9i2

    gg     = 1.0e0_rt + ee + ff
    dgg    = dee + dff

    term    = (term1 + term2)/gg
    dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


    ! the rates
    rev      = 1.81_rt * exp(-18.572_rt*tf%t9i)
    drevdt   = rev*18.572_rt*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    !drrdd = term

  end subroutine rate_mg24ap



  subroutine rate_al27pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                        dd,ddd,ee,dee,ff,dff,gg,dgg

    !$gpu

    ! al27(p,g)si28
    ! champagne 1996

    aa  = 1.32e09_rt * tf%t9i23 * exp(-23.26_rt*tf%t9i13)
    daa = aa*(-twoth*tf%t9i + oneth*23.26_rt*tf%t9i43)

    bb  = 3.22e-10_rt * tf%t9i32 * exp(-0.836_rt*tf%t9i)*0.17_rt
    dbb = bb*(-1.5e0_rt*tf%t9i + 0.836_rt*tf%t9i2)

    cc  = 1.74e00_rt * tf%t9i32 * exp(-2.269_rt*tf%t9i)
    dcc = cc*(-1.5e0_rt*tf%t9i + 2.269_rt*tf%t9i2)

    dd  = 9.92e00_rt * tf%t9i32 * exp(-2.492_rt*tf%t9i)
    ddd = dd*(-1.5e0_rt*tf%t9i + 2.492_rt*tf%t9i2)

    ee  = 4.29e01_rt * tf%t9i32 * exp(-3.273_rt*tf%t9i)
    dee = ee*(-1.5e0_rt*tf%t9i + 3.273_rt*tf%t9i2)

    ff  = 1.34e02_rt * tf%t9i32 * exp(-3.654_rt*tf%t9i)
    dff = ff*(-1.5e0_rt*tf%t9i + 3.654_rt*tf%t9i2)

    gg  = 1.77e04_rt * (tf%t9**0.53_rt) * exp(-4.588_rt*tf%t9i)
    dgg = gg*(0.53_rt*tf%t9i + 4.588_rt*tf%t9i2)

    term    = aa + bb + cc + dd + ee + ff + gg
    dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.13e11_rt * tf%t932 * exp(-134.434_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 134.434_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_al27pg



  subroutine rate_al27pg_old(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                        ff,dff,gg,dgg,hh,dhh,xx,dxx,yy,dyy,zz,dzz,pp, &
                        rev,drevdt

    real(rt), parameter :: rc147 = 0.1e0_rt
    real(rt), parameter :: q1    = 1.0e0_rt/0.024025e0_rt

    !$gpu

    ! 27al(p,g)si28  cf88
    aa  = 1.67e08_rt * tf%t9i23 * exp(-23.261_rt*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*23.261_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb  = 1.0e0_rt + 0.018_rt*tf%t913 + 5.81_rt*tf%t923 + 0.728_rt*tf%t9 &
         + 27.31_rt*tf%t943 + 8.71_rt*tf%t953
    dbb = oneth*0.018_rt*tf%t9i23 + twoth*5.81_rt*tf%t9i13 + 0.728_rt &
         + fourth*27.31_rt*tf%t913 + fiveth*8.71_rt*tf%t923

    cc  = aa*bb
    dcc = daa*bb + aa*dbb

    dd  = 2.20e00_rt * tf%t9i32 * exp(-2.269_rt*tf%t9i)
    ddd = dd*(-1.5e0_rt*tf%t9i + 2.269_rt*tf%t9i2)

    ee  = 1.22e01_rt * tf%t9i32 * exp(-2.491_rt*tf%t9i)
    dee = ee*(-1.5e0_rt*tf%t9i + 2.491_rt*tf%t9i2)

    ff  =  1.50e04_rt * tf%t9 * exp(-4.112_rt*tf%t9i)
    dff = ff*(tf%t9i + 4.112_rt*tf%t9i2)

    gg  = rc147 * 6.50e-10_rt * tf%t9i32 * exp(-0.853_rt*tf%t9i)
    dgg = gg*(-1.5e0_rt*tf%t9i + 0.853_rt*tf%t9i2)

    hh  = rc147 * 1.63e-10_rt * tf%t9i32 * exp(-1.001_rt*tf%t9i)
    dhh = hh*(-1.5e0_rt*tf%t9i + 1.001_rt*tf%t9i2)

    xx     = oneth*exp(-9.792_rt*tf%t9i)
    dxx    = xx*9.792_rt*tf%t9i2

    yy     =  twoth * exp(-11.773_rt*tf%t9i)
    dyy    = yy*11.773_rt*tf%t9i2

    zz     = 1.0e0_rt + xx + yy
    dzz    = dxx + dyy

    pp      = 1.0e0_rt/zz
    term    = (cc + dd + ee + ff + gg + hh)*pp
    dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.13e11_rt*tf%t932*exp(-134.434_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 134.434_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_al27pg_old



  subroutine rate_si28ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! si28(a,g)s32
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 6.340e-2_rt*z + 2.541e-3_rt*z2 - 2.900e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0
    else
       daa   = 6.340e-2_rt + 2.0e0_rt*2.541e-3_rt*tf%t9 - 3.0e0_rt*2.900e-4_rt*tf%t92
    end if

    term    = 4.82e22_rt * tf%t9i23 * exp(-61.015_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 61.015_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))
  
    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 6.461e10_rt * tf%t932 * exp(-80.643_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 80.643_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_si28ag



  subroutine rate_si28ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! si28(a,p)p31
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 2.798e-3_rt*z + 2.763e-3_rt*z2 - 2.341e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 2.798e-3_rt + 2.0e0_rt*2.763e-3_rt*tf%t9 - 3.0e0_rt*2.341e-4_rt*tf%t92
    end if

    term    = 4.16e13_rt * tf%t9i23 * exp(-25.631_rt * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*25.631_rt*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 0.5825e0_rt * exp(-22.224_rt*tf%t9i)
    drevdt   = rev*22.224_rt*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    !drrdd = term

  end subroutine rate_si28ap



  subroutine rate_p31pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! p31(p,g)s32
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.928e-1_rt*z - 1.540e-2_rt*z2 + 6.444e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 1.928e-1_rt - 2.0e0_rt*1.540e-2_rt*tf%t9 + 3.0e0_rt*6.444e-4_rt*tf%t92
    end if

    term    = 1.08e16_rt * tf%t9i23 * exp(-27.042_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 27.042_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 3.764e10_rt * tf%t932 * exp(-102.865_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 102.865_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_p31pg



  subroutine rate_s32ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! s32(a,g)ar36
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 4.913e-2_rt*z + 4.637e-3_rt*z2 - 4.067e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 4.913e-2_rt + 2.0e0_rt*4.637e-3_rt*tf%t9 - 3.0e0_rt*4.067e-4_rt*tf%t92
    end if

    term    = 1.16e24_rt * tf%t9i23 * exp(-66.690_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 66.690_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 6.616e10_rt * tf%t932 * exp(-77.080_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 77.080_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_s32ag



  subroutine rate_s32ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! s32(a,p)cl35
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.041e-1_rt*z - 1.368e-2_rt*z2 + 6.969e-4_rt*z3
    if (z .eq. 10) then
       daa = 0.0e0_rt
    else
       daa   = 1.041e-1_rt - 2.0e0_rt*1.368e-2_rt*tf%t9 + 3.0e0_rt*6.969e-4_rt*tf%t92
    end if

    term    = 1.27e16_rt * tf%t9i23 * exp(-31.044_rt * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*31.044_rt*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 1.144_rt * exp(-21.643_rt*tf%t9i)
    drevdt   = rev*21.643_rt*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    !drrdd = term

  end subroutine rate_s32ap



  subroutine rate_cl35pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt

    !$gpu

    ! cl35(p,g)ar36
    aa    = 1.0e0_rt + 1.761e-1_rt*tf%t9 - 1.322e-2_rt*tf%t92 + 5.245e-4_rt*tf%t93
    daa   = 1.761e-1_rt - 2.0e0_rt*1.322e-2_rt*tf%t9 + 3.0e0_rt*5.245e-4_rt*tf%t92
  

    term    =  4.48e16_rt * tf%t9i23 * exp(-29.483_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 29.483_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 7.568e10_rt*tf%t932*exp(-98.722_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 98.722_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_cl35pg


  subroutine rate_ar36ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ar36(a,g)ca40
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.458e-1_rt*z - 1.069e-2_rt*z2 + 3.790e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 1.458e-1_rt - 2.0e0_rt*1.069e-2_rt*tf%t9 + 3.0e0_rt*3.790e-4_rt*tf%t92
    end if

    term    = 2.81e30_rt * tf%t9i23 * exp(-78.271_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 78.271_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 6.740e10_rt * tf%t932 * exp(-81.711_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 81.711_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_ar36ag



  subroutine rate_ar36ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ar36(a,p)k39
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 4.826e-3_rt*z - 5.534e-3_rt*z2 + 4.021e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 4.826e-3_rt - 2.0e0_rt*5.534e-3_rt*tf%t9 + 3.0e0_rt*4.021e-4_rt*tf%t92
    end if

    term    = 2.76e13_rt * tf%t9i23 * exp(-34.922_rt * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*34.922_rt*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 1.128_rt*exp(-14.959_rt*tf%t9i)
    drevdt   = rev*14.959_rt*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    !drrdd = term

  end subroutine rate_ar36ap



  subroutine rate_k39pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! k39(p,g)ca40
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.622e-1_rt*z - 1.119e-2_rt*z2 + 3.910e-4_rt*z3
    if (z .eq. 10) then
       daa = 0.0e0_rt
    else
       daa   = 1.622e-1_rt - 2.0e0_rt*1.119e-2_rt*tf%t9 + 3.0e0_rt*3.910e-4_rt*tf%t92
    end if

    term    = 4.09e16_rt * tf%t9i23 * exp(-31.727_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 31.727_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 7.600e10_rt * tf%t932 * exp(-96.657_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 96.657_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_k39pg



  subroutine rate_ca40ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ca40(a,g)ti44
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.650e-2_rt*z + 5.973e-3_rt*z2 - 3.889e-04_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 1.650e-2_rt + 2.0e0_rt*5.973e-3_rt*tf%t9 - 3.0e0_rt*3.889e-4_rt*tf%t92
    end if

    term    = 4.66e24_rt * tf%t9i23 * exp(-76.435_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 76.435_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 6.843e10_rt * tf%t932 * exp(-59.510_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 59.510_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_ca40ag



  subroutine rate_ca40ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ca40(a,p)sc43
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt - 1.206e-2_rt*z + 7.753e-3_rt*z2 - 5.071e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = -1.206e-2_rt + 2.0e0_rt*7.753e-3_rt*tf%t9 - 3.0e0_rt*5.071e-4_rt*tf%t92
    end if

    term    = 4.54e14_rt * tf%t9i23 * exp(-32.177_rt * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*32.177_rt*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 2.229_rt * exp(-40.966_rt*tf%t9i)
    drevdt   = rev*40.966_rt*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    !drrdd = term

  end subroutine rate_ca40ap



  subroutine rate_sc43pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! sc43(p,g)ca40
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.023e-1_rt*z - 2.242e-3_rt*z2 - 5.463e-5_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 1.023e-1_rt - 2.0e0_rt*2.242e-3_rt*tf%t9 - 3.0e0_rt*5.463e-5_rt*tf%t92
    end if

    term    = 3.85e16_rt * tf%t9i23 * exp(-33.234_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 33.234_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.525e11_rt * tf%t932 * exp(-100.475_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 100.475_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_sc43pg



  subroutine rate_ti44ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ti44(a,g)cr48
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.066e-1_rt*z - 1.102e-2_rt*z2 + 5.324e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 1.066e-1_rt - 2.0e0_rt*1.102e-2_rt*tf%t9 + 3.0e0_rt*5.324e-4_rt*tf%t92
    end if

    term    = 1.37e26_rt * tf%t9i23 * exp(-81.227_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 81.227_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 6.928e10_rt*tf%t932*exp(-89.289_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 89.289_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_ti44ag



  subroutine rate_ti44ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! ti44(a,p)v47
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 2.655e-2_rt*z - 3.947e-3_rt*z2 + 2.522e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 2.655e-2_rt - 2.0e0_rt*3.947e-3_rt*tf%t9 + 3.0e0_rt*2.522e-4_rt*tf%t92
    end if

    term    = 6.54e20_rt * tf%t9i23 * exp(-66.678_rt * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*66.678_rt*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    rev      = 1.104_rt * exp(-4.723_rt*tf%t9i)
    drevdt   = rev*4.723_rt*tf%t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !dfrdd = rev * term

    rr    = den * term
    drrdt = den * dtermdt * 1.0e-9_rt
    !drrdd = term

  end subroutine rate_ti44ap



  subroutine rate_v47pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! v47(p,g)cr48
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 9.979e-2_rt*z - 2.269e-3_rt*z2 - 6.662e-5_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 9.979e-2_rt - 2.0e0_rt*2.269e-3_rt*tf%t9 - 3.0e0_rt*6.662e-5_rt*tf%t92
    end if

    term    = 2.05e17_rt * tf%t9i23 * exp(-35.568_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 35.568_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 7.649e10_rt*tf%t932*exp(-93.999_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 93.999_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_v47pg



  subroutine rate_cr48ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! cr48(a,g)fe52
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 6.325e-2_rt*z - 5.671e-3_rt*z2 + 2.848e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 6.325e-2_rt - 2.0e0_rt*5.671e-3_rt*tf%t9 + 3.0e0_rt*2.848e-4_rt*tf%t92
    end if

    term    = 1.04e23_rt * tf%t9i23 * exp(-81.420_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 81.420_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 7.001e10_rt * tf%t932 * exp(-92.177_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 92.177_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_cr48ag



  subroutine rate_cr48ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! cr48(a,p)mn51
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.384e-2_rt*z + 1.081e-3_rt*z2 - 5.933e-5_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 1.384e-2_rt + 2.0e0_rt*1.081e-3_rt*tf%t9 - 3.0e0_rt*5.933e-5_rt*tf%t92
    end if

    term    = 1.83e26_rt * tf%t9i23 * exp(-86.741_rt * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*86.741_rt*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 0.6087_rt*exp(-6.510_rt*tf%t9i)
    drevdt   = rev*6.510_rt*tf%t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = rev * term

  end subroutine rate_cr48ap



  subroutine rate_mn51pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! mn51(p,g)fe52
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 8.922e-2_rt*z - 1.256e-3_rt*z2 - 9.453e-5_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 8.922e-2_rt - 2.0e0_rt*1.256e-3_rt*tf%t9 - 3.0e0_rt*9.453e-5_rt*tf%t92
    end if

    term    = 3.77e17_rt * tf%t9i23 * exp(-37.516_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 37.516_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.150e11_rt*tf%t932*exp(-85.667_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 85.667_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_mn51pg



  subroutine rate_fe52ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! fe52(a,g)ni56
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 7.846e-2_rt*z - 7.430e-3_rt*z2 + 3.723e-4_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 7.846e-2_rt - 2.0e0_rt*7.430e-3_rt*tf%t9 + 3.0e0_rt*3.723e-4_rt*tf%t92
    end if

    term    = 1.05e27_rt * tf%t9i23 * exp(-91.674_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 91.674_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 7.064e10_rt*tf%t932*exp(-92.850_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 92.850_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_fe52ag



  subroutine rate_fe52ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! fe52(a,p)co55
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 1.367e-2_rt*z + 7.428e-4_rt*z2 - 3.050e-5_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 1.367e-2_rt + 2.0e0_rt*7.428e-4_rt*tf%t9 - 3.0e0_rt*3.050e-5_rt*tf%t92
    end if

    term    = 1.30e27_rt * tf%t9i23 * exp(-91.674_rt * tf%t9i13 * aa)
    dtermdt = -twoth*term*tf%t9i + term*91.674_rt*tf%t9i13*(oneth*tf%t9i*aa - daa)


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 0.4597_rt*exp(-9.470_rt*tf%t9i)
    drevdt   = rev*9.470_rt*tf%t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = rev * term

  end subroutine rate_fe52ap



  subroutine rate_co55pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,rev,drevdt,z,z2,z3

    !$gpu

    ! co55(p,g)ni56
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 9.894e-2_rt*z - 3.131e-3_rt*z2 - 2.160e-5_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 9.894e-2_rt - 2.0e0_rt*3.131e-3_rt*tf%t9 - 3.0e0_rt*2.160e-5_rt*tf%t92
    end if

    term    = 1.21e18_rt * tf%t9i23 * exp(-39.604_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 39.604_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.537e11_rt*tf%t932*exp(-83.382_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 83.382_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_co55pg



  subroutine rate_pp(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,aa,daa,bb,dbb

    !$gpu

    ! p(p,e+nu)d
    if (tf%t9 .le. 3.0_rt) then
       aa   = 4.01e-15_rt * tf%t9i23 * exp(-3.380e0_rt*tf%t9i13)
       daa  = aa*(-twoth*tf%t9i + oneth*3.380e0_rt*tf%t9i43)

       bb   = 1.0e0_rt + 0.123e0_rt*tf%t913 + 1.09e0_rt*tf%t923 + 0.938e0_rt*tf%t9
       dbb  = oneth*0.123e0_rt*tf%t9i23 + twoth*1.09e0_rt*tf%t9i13 + 0.938e0_rt

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

    else
       term    = 1.1581136e-15_rt
       dtermdt = 0.0e0_rt
    end if

    ! rate
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rr    = 0.0e0_rt
    drrdt = 0.0e0_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_pp



  subroutine rate_png(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa

    !$gpu

    ! p(n,g)d
    ! smith,kawano,malany 1992

    aa      = 1.0e0_rt - 0.8504_rt*tf%t912 + 0.4895_rt*tf%t9 &
         - 0.09623_rt*tf%t932 + 8.471e-3_rt*tf%t92 &
         - 2.80e-4_rt*tf%t952

    daa     =  -0.5e0_rt*0.8504_rt*tf%t9i12 + 0.4895_rt &
         - 1.5e0_rt*0.09623_rt*tf%t912 + 2.0e0_rt*8.471e-3_rt*tf%t9 &
         - 2.5e0_rt*2.80e-4_rt*tf%t932

    term    = 4.742e4_rt * aa
    dtermdt = 4.742e4_rt * daa


    ! wagoner,schramm 1977
    !      aa      = 1.0e0_rt - 0.86_rt*tf%t912 + 0.429_rt*tf%t9
    !      daa     =  -0.5e0_rt*0.86_rt*tf%t9i12 + 0.429

    !      term    = 4.4e4_rt * aa
    !      dtermdt = 4.4e4_rt * daa



    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 4.71e09_rt * tf%t932 * exp(-25.82_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 25.82_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_png



  subroutine rate_dpg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb

    !$gpu

    ! d(p,g)he3
    aa      = 2.24e03_rt * tf%t9i23 * exp(-3.720_rt*tf%t9i13)
    daa     = aa*(-twoth*tf%t9i + oneth*3.720_rt*tf%t9i43)

    bb      = 1.0e0_rt + 0.112_rt*tf%t913 + 3.38_rt*tf%t923 + 2.65_rt*tf%t9
    dbb     = oneth*0.112_rt*tf%t9i23 + twoth*3.38_rt*tf%t9i13 + 2.65_rt

    term    = aa * bb
    dtermdt = daa * bb + aa * dbb


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.63e10_rt * tf%t932 * exp(-63.750_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 63.750_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_dpg



  subroutine rate_he3ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt

    !$gpu

    ! he3(n,g)he4
    term    = 6.62_rt * (1.0e0_rt + 905.0_rt*tf%t9)
    dtermdt = 5.9911e3_rt

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 2.61e10_rt * tf%t932 * exp(-238.81_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 238.81_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_he3ng



  subroutine rate_he3he3(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb

    !$gpu

    ! he3(he3,2p)he4
    aa   = 6.04e10_rt * tf%t9i23 * exp(-12.276_rt*tf%t9i13)
    daa  = aa*(-twoth*tf%t9i + oneth*12.276_rt*tf%t9i43)

    bb   = 1.0e0_rt + 0.034_rt*tf%t913 - 0.522_rt*tf%t923 - 0.124_rt*tf%t9 &
         + 0.353_rt*tf%t943 + 0.213_rt*tf%t953
    dbb  = oneth*0.034_rt*tf%t9i23 - twoth*0.522_rt*tf%t9i13 - 0.124_rt &
         + fourth*0.353_rt*tf%t913 + fiveth*0.213_rt*tf%t923

    term    = aa * bb
    dtermdt = daa*bb + aa*dbb

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 3.39e-10_rt * tf%t9i32 * exp(-149.230_rt*tf%t9i)
    drevdt   = rev*(-1.5e0_rt*tf%t9i + 149.230_rt*tf%t9i2)

    rr    = den * den * rev * term
    drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 2.0e0_rt * den * rev * term

  end subroutine rate_he3he3



  subroutine rate_he3he4(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,t9a,dt9a, &
                        t9a13,dt9a13,t9a56,dt9a56,zz

    !$gpu

    ! he3(he4,g)be7
    aa      = 1.0e0_rt + 0.0495_rt*tf%t9
    daa     = 0.0495_rt

    zz      = 1.0e0_rt/aa
    t9a     = tf%t9*zz
    dt9a    = (1.0e0_rt - t9a*daa)*zz

    zz      = dt9a/t9a
    t9a13   = t9a**oneth
    dt9a13  = oneth*t9a13*zz

    t9a56   = t9a**fivsix
    dt9a56  = fivsix*t9a56*zz

    term    = 5.61e6_rt * t9a56 * tf%t9i32 * exp(-12.826_rt/t9a13)
    dtermdt = term*(dt9a56/t9a56 - 1.5e0_rt*tf%t9i &
         + 12.826_rt/t9a13**2 * dt9a13)

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.11e10_rt * tf%t932 * exp(-18.423_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 18.423_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_he3he4



  subroutine rate_c12pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee

    real(rt), parameter :: q1 = 1.0e0_rt/2.25e0_rt

    !$gpu

    ! c12(p,g)13n
    aa   = 2.04e07_rt * tf%t9i23 * exp(-13.69_rt*tf%t9i13 - tf%t92*q1)
    daa  = aa*(-twoth*tf%t9i + oneth*13.69_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb   = 1.0e0_rt + 0.03_rt*tf%t913 + 1.19_rt*tf%t923 + 0.254_rt*tf%t9 &
         + 2.06_rt*tf%t943 + 1.12_rt*tf%t953
    dbb  = oneth*0.03_rt*tf%t9i23 + twoth*1.19_rt*tf%t9i13 + 0.254_rt &
         + fourth*2.06_rt*tf%t913 + fiveth*1.12_rt*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 1.08e05_rt * tf%t9i32 * exp(-4.925_rt*tf%t9i)
    ddd  = dd*(-1.5e0_rt*tf%t9i + 4.925_rt*tf%t9i2)

    ee   = 2.15e05_rt * tf%t9i32 * exp(-18.179_rt*tf%t9i)
    dee  = ee*(-1.5e0_rt*tf%t9i + 18.179_rt*tf%t9i2)

    term    = cc + dd + ee
    dtermdt = dcc + ddd + dee

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    dfrdd = term

    rev      = 8.84e09_rt * tf%t932 * exp(-22.553_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 22.553_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    drrdd = 0.0e0_rt

  end subroutine rate_c12pg



  subroutine rate_n14pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee

    real(rt), parameter :: q1 = 1.0e0_rt/10.850436e0_rt

    !$gpu

    ! n14(p,g)o15
    aa  = 4.90e07_rt * tf%t9i23 * exp(-15.228_rt*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*15.228_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb   = 1.0e0_rt + 0.027_rt*tf%t913 - 0.778_rt*tf%t923 - 0.149_rt*tf%t9 &
         + 0.261_rt*tf%t943 + 0.127_rt*tf%t953
    dbb  = oneth*0.027_rt*tf%t9i23 - twoth*0.778_rt*tf%t9i13 - 0.149_rt &
         + fourth*0.261_rt*tf%t913 + fiveth*0.127_rt*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 2.37e03_rt * tf%t9i32 * exp(-3.011_rt*tf%t9i)
    ddd  = dd*(-1.5e0_rt*tf%t9i + 3.011_rt*tf%t9i2)

    ee   = 2.19e04_rt * exp(-12.530_rt*tf%t9i)
    dee  = ee*12.530_rt*tf%t9i2

    term    = cc + dd + ee
    dtermdt = dcc + ddd + dee

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev    = 2.70e10_rt * tf%t932 * exp(-84.678_rt*tf%t9i)
    drevdt = rev*(1.5e0_rt*tf%t9i + 84.678_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_n14pg



  subroutine rate_n15pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,ff,dff

    real(rt), parameter :: q1 = 1.0e0_rt/0.2025e0_rt

    !$gpu

    ! n15(p,g)o16
    aa  = 9.78e08_rt * tf%t9i23 * exp(-15.251_rt*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*15.251_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb   = 1.0e0_rt  + 0.027_rt*tf%t913 + 0.219_rt*tf%t923 + 0.042_rt*tf%t9 &
         + 6.83_rt*tf%t943 + 3.32_rt*tf%t953
    dbb  = oneth*0.027_rt*tf%t9i23 + twoth*0.219_rt*tf%t9i13 + 0.042_rt &
         + fourth*6.83_rt*tf%t913 + fiveth*3.32_rt*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 1.11e04_rt*tf%t9i32*exp(-3.328_rt*tf%t9i)
    ddd  = dd*(-1.5e0_rt*tf%t9i + 3.328_rt*tf%t9i2)

    ee   = 1.49e04_rt*tf%t9i32*exp(-4.665_rt*tf%t9i)
    dee  = ee*(-1.5e0_rt*tf%t9i + 4.665_rt*tf%t9i2)

    ff   = 3.8e06_rt*tf%t9i32*exp(-11.048_rt*tf%t9i)
    dff  = ff*(-1.5e0_rt*tf%t9i + 11.048_rt*tf%t9i2)

    term    = cc + dd + ee + ff
    dtermdt = dcc + ddd + dee + dff

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 3.62e10_rt * tf%t932 * exp(-140.734_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 140.734_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_n15pg



  subroutine rate_n15pa(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg

    real(rt), parameter :: theta = 0.1e0_rt
    real(rt), parameter :: q1    = 1.0e0_rt/0.272484e0_rt

    !$gpu

    ! n15(p,a)c12
    aa  = 1.08e12_rt*tf%t9i23*exp(-15.251_rt*tf%t9i13 - tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*15.251_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb   = 1.0e0_rt + 0.027_rt*tf%t913 + 2.62_rt*tf%t923 + 0.501_rt*tf%t9 &
         + 5.36_rt*tf%t943 + 2.60_rt*tf%t953
    dbb  = oneth*0.027_rt*tf%t9i23 + twoth*2.62_rt*tf%t9i13 + 0.501_rt &
         + fourth*5.36_rt*tf%t913 + fiveth*2.60_rt*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 1.19e08_rt * tf%t9i32 * exp(-3.676_rt*tf%t9i)
    ddd  = dd*(-1.5e0_rt*tf%t9i + 3.676_rt*tf%t9i2)

    ee   = 5.41e08_rt * tf%t9i12 * exp(-8.926_rt*tf%t9i)
    dee  = ee*(-0.5e0_rt*tf%t9i + 8.926_rt*tf%t9i2)

    ff   = theta * 4.72e08_rt * tf%t9i32 * exp(-7.721_rt*tf%t9i)
    dff  = ff*(-1.5e0_rt*tf%t9i + 7.721_rt*tf%t9i2)

    gg   = theta * 2.20e09_rt * tf%t9i32 * exp(-11.418_rt*tf%t9i)
    dgg  = gg*(-1.5e0_rt*tf%t9i + 11.418_rt*tf%t9i2)

    term    = cc + dd + ee + ff + gg
    dtermdt = dcc + ddd + dee + dff + dgg

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 7.06e-01_rt*exp(-57.625_rt*tf%t9i)
    drevdt   = rev*57.625_rt*tf%t9i2

    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = rev * term

  end subroutine rate_n15pa



  subroutine rate_o16pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,zz

    !$gpu

    ! o16(p,g)f17
    aa  = exp(-0.728_rt*tf%t923)
    daa = -twoth*aa*0.728_rt*tf%t9i13

    bb  = 1.0e0_rt + 2.13_rt * (1.0e0_rt - aa)
    dbb = -2.13_rt*daa

    cc  = tf%t923 * bb
    dcc = twoth*cc*tf%t9i + tf%t923*dbb

    dd   = exp(-16.692_rt*tf%t9i13)
    ddd  = oneth*dd*16.692_rt*tf%t9i43

    zz   = 1.0e0_rt/cc
    ee   = dd*zz
    dee  = (ddd - ee*dcc)*zz

    term    = 1.50e08_rt * ee
    dtermdt = 1.50e08_rt * dee


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 3.03e09_rt*tf%t932*exp(-6.968_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 6.968_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_o16pg



  subroutine rate_n14ag(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                        cc,dcc,dd,ddd,ee,dee,ff,dff

    real(rt), parameter :: q1 = 1.0e0_rt/0.776161e0_rt

    !$gpu

    ! n14(a,g)f18
    aa  = 7.78e09_rt * tf%t9i23 * exp(-36.031_rt*tf%t9i13- tf%t92*q1)
    daa = aa*(-twoth*tf%t9i + oneth*36.031_rt*tf%t9i43 - 2.0e0_rt*tf%t9*q1)

    bb   = 1.0e0_rt + 0.012_rt*tf%t913 + 1.45_rt*tf%t923 + 0.117_rt*tf%t9 &
         + 1.97_rt*tf%t943 + 0.406_rt*tf%t953
    dbb  = oneth*0.012_rt*tf%t9i23 + twoth*1.45_rt*tf%t9i13 + 0.117_rt &
         + fourth*1.97_rt*tf%t913 + fiveth*0.406_rt*tf%t923

    cc   = aa * bb
    dcc  = daa*bb + aa*dbb

    dd   = 2.36e-10_rt * tf%t9i32 * exp(-2.798_rt*tf%t9i)
    ddd  = dd*(-1.5e0_rt*tf%t9i + 2.798_rt*tf%t9i2)

    ee   = 2.03_rt * tf%t9i32 * exp(-5.054_rt*tf%t9i)
    dee  = ee*(-1.5e0_rt*tf%t9i + 5.054_rt*tf%t9i2)

    ff   = 1.15e04_rt * tf%t9i23 * exp(-12.310_rt*tf%t9i)
    dff  = ff*(-twoth*tf%t9i + 12.310_rt*tf%t9i2)

    term    = cc + dd + ee + ff
    dtermdt = dcc + ddd + dee + dff

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 5.42e10_rt * tf%t932 * exp(-51.236_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 51.236_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_n14ag



  subroutine rate_fe52ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,tq2

    !$gpu

    ! fe52(n,g)fe53
    tq2     = tf%t9 - 0.348e0_rt
    term    = 9.604e05_rt * exp(-0.0626_rt*tq2)
    dtermdt = -term*0.0626_rt

    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 2.43e09_rt * tf%t932 * exp(-123.951_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 123.951_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_fe52ng



  subroutine rate_fe53ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,tq1,tq10,dtq10,tq2

    !$gpu

    ! fe53(n,g)fe54
    tq1   = tf%t9/0.348_rt
    tq10  = tq1**0.10_rt
    dtq10 = 0.1e0_rt*tq10/(0.348_rt*tq1)
    tq2   = tf%t9 - 0.348e0_rt

    term    = 1.817e06_rt * tq10 * exp(-0.06319_rt*tq2)
    dtermdt = term/tq10*dtq10 - term*0.06319_rt

    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 1.56e11_rt * tf%t932 * exp(-155.284_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 155.284_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_fe53ng



  subroutine rate_fe54ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den, fr, dfrdt, dfrdd, rr, drrdt, drrdd
    type (tf_t)      :: tf

    real(rt) :: aa, daa, bb, dbb, term, dtermdt

    !$gpu

    ! fe54(n,g)fe55
    aa   =  2.307390e01_rt - 7.931795e-02_rt * tf%t9i + 7.535681e00_rt * tf%t9i13 &
         - 1.595025e01_rt * tf%t913 + 1.377715e00_rt * tf%t9 - 1.291479e-01_rt * tf%t953 &
         + 6.707473e00_rt * log(tf%t9)

    daa  =  7.931795e-02_rt * tf%t9i2 - oneth * 7.535681e00_rt * tf%t9i43 &
         - oneth * 1.595025e01_rt *tf%t9i23 + 1.377715e00_rt - fiveth * 1.291479e-01_rt *tf%t923 &
         + 6.707473e00_rt * tf%t9i

    if (aa .lt. 200.0_rt) then
       term    = exp(aa)
       dtermdt = term*daa*1.0e-9_rt
    else
       term    = exp(200.0e0_rt)
       dtermdt = 0.0e0_rt
    end if

    bb  = 4.800293e09_rt * tf%t932 * exp(-1.078986e02_rt * tf%t9i)
    dbb = bb*(1.5e0_rt*tf%t9i + 1.078986e02_rt * tf%t9i2)

    ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0e-9_rt
    !drrdd = 0.0e0_rt

    ! forward rate
    !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den

  end subroutine rate_fe54ng



  subroutine rate_fe54pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: term,dtermdt,rev,drevdt,aa,daa,z,z2,z3

    !$gpu

    ! fe54(p,g)co55
    z     = min(tf%t9,10.0e0_rt)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0e0_rt + 9.593e-2_rt*z - 3.445e-3_rt*z2 + 8.594e-5_rt*z3
    if (z .eq. 10.0_rt) then
       daa = 0.0e0_rt
    else
       daa   = 9.593e-2_rt - 2.0e0_rt*3.445e-3_rt*tf%t9 + 3.0e0_rt*8.594e-5_rt*tf%t92
    end if

    term    = 4.51e17_rt * tf%t9i23 * exp(-38.483_rt * tf%t9i13 * aa)
    dtermdt = term*(-twoth*tf%t9i + 38.483_rt*tf%t9i13*(oneth*tf%t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0e-9_rt
    !dfrdd = term

    rev      = 2.400e09_rt * tf%t932 * exp(-58.605_rt*tf%t9i)
    drevdt   = rev*(1.5e0_rt*tf%t9i + 58.605_rt*tf%t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0e-9_rt
    !drrdd = 0.0e0_rt

  end subroutine rate_fe54pg




  subroutine rate_fe54ap(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: aa,daa,bb,dbb,term,dtermdt

    !$gpu

    ! fe54(a,p)co57
    aa   =  3.97474900e01_rt - 6.06543100e00_rt * tf%t9i + 1.63239600e02_rt * tf%t9i13 &
         - 2.20457700e02_rt * tf%t913 + 8.63980400e00_rt * tf%t9 - 3.45841300e-01_rt * tf%t953 &
         + 1.31464200e02_rt * log(tf%t9)

    daa  =  6.06543100e00_rt * tf%t9i2 - oneth * 1.63239600e02_rt * tf%t9i43 &
         - oneth * 2.20457700e02_rt * tf%t9i23 + 8.63980400e00_rt - fiveth * 3.45841300e-01_rt * tf%t923 &
         + 1.31464200e02_rt  * tf%t9i

    if (aa .lt. 200.0_rt) then
       term    = exp(aa)
       dtermdt = term*daa*1.0e-9_rt
    else
       term    = exp(200.0e0_rt)
       dtermdt = 0.0e0_rt
    end if

    bb  = 2.16896000e00_rt  * exp(-2.05631700e01_rt * tf%t9i)
    dbb = bb * 2.05631700e01_rt * tf%t9i2

    ! reverse rate
    !drrdd = term
    rr    = term*den
    drrdt = dtermdt*den

    ! forward rate
    fr    = rr*bb
    dfrdt = drrdt*bb + rr*dbb*1.0e-9_rt
    !dfrdd = drrdd*bb

  end subroutine rate_fe54ap



  subroutine rate_fe55ng(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: aa,daa,bb,dbb,term,dtermdt

    !$gpu

    ! fe55(n,g)fe56
    aa   =  1.954115e01_rt - 6.834029e-02_rt * tf%t9i + 5.379859e00_rt * tf%t9i13 &
         - 8.758150e00_rt * tf%t913 + 5.285107e-01_rt * tf%t9 - 4.973739e-02_rt  * tf%t953 &
         + 4.065564e00_rt  * log(tf%t9)

    daa  =  6.834029e-02_rt * tf%t9i2 - oneth * 5.379859e00_rt * tf%t9i43 &
         - oneth * 8.758150e00_rt * tf%t9i23 + 5.285107e-01_rt - fiveth * 4.973739e-02_rt  *tf%t923 &
         + 4.065564e00_rt  * tf%t9i

    if (aa .lt. 200.0_rt) then
       term    = exp(aa)
       dtermdt = term*daa*1.0e-9_rt
    else
       term    = exp(200.0e0_rt)
       dtermdt = 0.0e0_rt
    end if

    bb  = 7.684279e10_rt  * tf%t932 * exp(-1.299472e02_rt  * tf%t9i)
    dbb = bb*(1.5e0_rt*tf%t9i + 1.299472e02_rt * tf%t9i2)

    ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0e-9_rt
    !drrdd = 0.0e0_rt

    ! forward rate
    !dfrdd = term
    fr    = term*den
    dfrdt = dtermdt*den

  end subroutine rate_fe55ng




  subroutine rate_fe56pg(tf,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    implicit none

    real(rt) :: den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    type (tf_t)      :: tf

    real(rt) :: aa,daa,bb,dbb,term,dtermdt

    !$gpu

    ! fe56(p,g)co57

    aa   =  1.755960e02_rt - 7.018872e00_rt * tf%t9i + 2.800131e02_rt * tf%t9i13 &
         - 4.749343e02_rt * tf%t913 + 2.683860e01_rt * tf%t9 - 1.542324e00_rt  * tf%t953 &
         + 2.315911e02_rt  * log(tf%t9)

    daa  =  7.018872e00_rt * tf%t9i2 - oneth * 2.800131e02_rt * tf%t9i43 &
         - oneth * 4.749343e02_rt * tf%t9i23 + 2.683860e01_rt - fiveth * 1.542324e00_rt  *tf%t923 &
         + 2.315911e02_rt  * tf%t9i

    if (aa .lt. 200.0_rt) then
       term    = exp(aa)
       dtermdt = term*daa*1.0e-9_rt
    else
       term    = exp(200.0e0_rt)
       dtermdt = 0.0e0_rt
    end if

    bb  = 2.402486e09_rt * tf%t932 * exp(-6.995192e01_rt * tf%t9i)
    dbb = bb*(1.5e0_rt*tf%t9i + 6.995192e01_rt * tf%t9i2)


    ! reverse rate
    rr    = term*bb
    drrdt = dtermdt*bb + term*dbb*1.0e-9_rt
    !drrdd = 0.0e0_rt

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
    real(rt) :: btemp,bden,y56,ye,rn56ec,sn56ec

    real(rt) :: rnt(2),rne(2,14),t9,r,rfm,rf0, &
                        rf1,rf2,dfacm,dfac0,dfac1,dfac2, &
                        tfm,tf0,tf1,tf2,tfacm,tfac0,tfac1,tfac2

    !$gpu

    ! calculate ni56 electron capture and neutrino loss rates
    rn56ec = 0.0_rt
    sn56ec = 0.0_rt
    if ( (btemp .lt. 1.0e9_rt) .or. (bden*ye .lt. 1.0e6_rt)) return
    t9    = min(btemp,1.4e10_rt) * 1.0e-9_rt
    r     = max(6.0e0_rt,min(11.0e0_rt,log10(bden*ye)))
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
    rn56ec = 10.0e0_rt**rnt(1)
    sn56ec = 6.022548e23_rt * 1.60218e-6_rt * y56 * 10.0e0_rt**rnt(2)
   
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

    real(rt) :: etakep,temp,rpen,rnep,spenc,snepc

    integer          iflag
    real(rt) :: t9,t5,qn,etaef,etael,zetan,eta,etael2, &
                        etael3,etael4,f1l,f2l,f3l,f4l,f5l,f1g, &
                        f2g,f3g,f4g,f5g,exmeta,eta2,eta3,eta4, &
                        fac0,fac1,fac2,fac3,rie1,rie2,facv0,facv1, &
                        facv2,facv3,facv4,rjv1,rjv2,spen,snep, &
                        pi2,exeta,zetan2,f0,etael5,bktinv, &
                        qn1,ftinv,twoln,cmk5,cmk6,bk,pi,qn2,c2me, &
                        xmp,xmn,qndeca,tmean
    parameter        (qn1    = -2.0716446e-06_rt, &
         ftinv  = 1.0e0_rt/1083.9269e0_rt, &
         twoln  = 0.6931472e0_rt, &
         cmk5   = 1.3635675e-49_rt, &
         cmk6   = 2.2993864e-59_rt, &
         bk     = 1.38062e-16_rt, &
         pi     = 3.1415927e0_rt, &
         pi2    = pi * pi, &
         qn2    = 2.0716446e-06_rt, &
         c2me   = 8.1872665e-07_rt, &
         xmp    = 1.6726485e-24_rt, &
         xmn    = 1.6749543e-24_rt, &
         qndeca = 1.2533036e-06_rt, &
         tmean  = 886.7e0_rt)
    !     3                  tmean  = 935.14e0_rt)

    real(rt) third,sixth
    parameter        (third = 1.0e0_rt/3.0e0_rt, &
         sixth = 1.0e0_rt/6.0e0_rt)

    !$gpu

    ! tmean and qndeca are the mean lifetime and decay energy of the neutron
    ! xmp,xnp are masses of the p and n in grams.
    ! c2me is the constant used to convert the neutrino energy
    ! loss rate from mec2/s (as in the paper) to ergs/particle/sec.

    ! initialize
    rpen   = 0.0e0_rt
    rnep   = 0.0e0_rt
    spen   = 0.0e0_rt
    snep   = 0.0e0_rt
    t9     = temp * 1.0e-9_rt
    bktinv = 1.0e0_rt/(bk *temp)
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
    if (eta .le. 6.8e02_rt) then
       exeta = exp(eta)
    else
       exeta = 0.0e0_rt
    end if
    etael2 = etael*etael
    etael3 = etael2*etael
    etael4 = etael3*etael
    etael5 = etael4*etael
    zetan2 = zetan*zetan
    if (eta .le. 6.8e02_rt) then
       f0 = log(1.0e0_rt + exeta)
    else
       f0 = eta
    end if

    ! if eta le. 0._rt, the following fermi integrals apply
    f1l = exeta
    f2l = 2.0e0_rt   * f1l
    f3l = 6.0e0_rt   * f1l
    f4l = 24.0e0_rt  * f1l
    f5l = 120.0e0_rt * f1l

    ! if eta gt. 0._rt, the following fermi integrals apply:
    f1g = 0.0e0_rt
    f2g = 0.0e0_rt
    f3g = 0.0e0_rt
    f4g = 0.0e0_rt
    f5g = 0.0e0_rt
    if (eta .gt. 0.0_rt) then
       exmeta = dexp(-eta)
       eta2   = eta*eta
       eta3   = eta2*eta
       eta4   = eta3*eta
       f1g = 0.5e0_rt*eta2 + 2.0e0_rt - exmeta
       f2g = eta3*third + 4.0e0_rt*eta + 2.0e0_rt*exmeta
       f3g = 0.25e0_rt*eta4 + 0.5e0_rt*pi2*eta2 + 12.0e0_rt - 6.0e0_rt*exmeta
       f4g = 0.2e0_rt*eta4*eta + 2.0e0_rt*pi2*third*eta3 + 48.0e0_rt*eta &
            + 24.0e0_rt*exmeta
       f5g = eta4*eta2*sixth + 5.0e0_rt*sixth*pi2*eta4 &
            + 7.0e0_rt*sixth*pi2*eta2  + 240.0e0_rt -120.e0_rt*exmeta
    end if

    ! factors which are multiplied by the fermi integrals
    fac3 = 2.0e0_rt*zetan + 4.0e0_rt*etael
    fac2 = 6.0e0_rt*etael2 + 6.0e0_rt*etael*zetan + zetan2
    fac1 = 4.0e0_rt*etael3 + 6.0e0_rt*etael2*zetan + 2.0e0_rt*etael*zetan2
    fac0 = etael4 + 2.0e0_rt*zetan*etael3 + etael2*zetan2

    ! electron capture rates onto protons with no blocking
    rie1 = f4l + fac3*f3l + fac2*f2l + fac1*f1l + fac0*f0
    rie2 = f4g + fac3*f3g + fac2*f2g + fac1*f1g + fac0*f0

    ! neutrino emission rate for electron capture:
    facv4 = 5.0e0_rt*etael + 3.0e0_rt*zetan
    facv3 = 10.0e0_rt*etael2 + 12.0e0_rt*etael*zetan + 3.0e0_rt*zetan2
    facv2 = 10.0e0_rt*etael3 + 18.0e0_rt*etael2*zetan &
         + 9.0e0_rt*etael*zetan2 + zetan2*zetan
    facv1 = 5.0e0_rt*etael4 + 12.0e0_rt*etael3*zetan &
         + 9.0e0_rt*etael2*zetan2 + 2.0e0_rt*etael*zetan2*zetan
    facv0 = etael5 + 3.0e0_rt*etael4*zetan &
         + 3.0e0_rt*etael3*zetan2 + etael2*zetan2*zetan
    rjv1  = f5l + facv4*f4l + facv3*f3l &
         + facv2*f2l + facv1*f1l + facv0*f0
    rjv2  = f5g + facv4*f4g + facv3*f3g &
         + facv2*f2g + facv1*f1g + facv0*f0

    ! for electrons capture onto protons
    if (iflag.eq.2) go to 503
    if (eta.gt.0._rt) go to 505
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
503 if (eta.gt.0._rt) go to 507
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
