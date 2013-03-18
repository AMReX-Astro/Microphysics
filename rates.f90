module rates_module

  ! these rate routine come directly from the orignal public_aprox13.f90
  ! file.  Only those used in this network are kept.  We modify the calling
  ! sequence to take the tfactors as an argument, and remove the include
  ! implno statement.

  implicit none

contains

  subroutine rate_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                     dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2, &
                     zz,q1
    parameter        (q1 = 1.0d0/12.222016d0)


    ! c12(a,g)o16
    aa   = 1.0d0 + 0.0489d0*t9i23
    daa  = -twoth*0.0489d0*t9i53

    bb   = t92*aa*aa
    dbb  = 2.0d0*(bb*t9i + t92*aa*daa)

    cc   = exp(-32.120d0*t9i13 - t92*q1)
    dcc  = cc * (oneth*32.120d0*t9i43 - 2.0d0*t9*q1)

    dd   = 1.0d0 + 0.2654d0*t9i23
    ddd  = -twoth*0.2654d0*t9i53

    ee   = t92*dd*dd
    dee  = 2.0d0*(ee*t9i + t92*dd*ddd)

    ff   = exp(-32.120d0*t9i13)
    dff  = ff * oneth*32.120d0*t9i43

    gg   = 1.25d3 * t9i32 * exp(-27.499*t9i)
    dgg  = gg*(-1.5d0*t9i + 27.499*t9i2)

    hh   = 1.43d-2 * t95 * exp(-15.541*t9i)
    dhh  = hh*(5.0d0*t9i + 15.541*t9i2)

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
    dfrdd = term

    rev    = 5.13d10 * t932 * exp(-83.111*t9i)
    drevdt = rev*(1.5d0*t9i + 83.111*t9i2)

    rr     = rev * term
    drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9
    drrdd  = 0.0d0

    return
  end subroutine rate_c12ag


  subroutine rate_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac, &
                     drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                     ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1,rc28, &
                     q1,q2
    parameter        (rc28   = 0.1d0, &
                      q1     = 1.0d0/0.009604d0, &
                      q2     = 1.0d0/0.055225d0)

    ! triple alfa to c12
    ! this is a(a,g)be8
    aa    = 7.40d+05 * t9i32 * exp(-1.0663*t9i)
    daa   = aa*(-1.5d0*t9i  + 1.0663*t9i2)
    
    bb    = 4.164d+09 * t9i23 * exp(-13.49*t9i13 - t92*q1)
    dbb   = bb*(-twoth*t9i + oneth*13.49*t9i43 - 2.0d0*t9*q1)
    
    cc    = 1.0d0 + 0.031*t913 + 8.009*t923 + 1.732*t9 &
         + 49.883*t943 + 27.426*t953
    dcc   = oneth*0.031*t9i23 + twoth*8.009*t9i13 + 1.732 &
         + fourth*49.883*t913 + fiveth*27.426*t923

    r2abe    = aa + bb * cc
    dr2abedt = daa + dbb*cc + bb*dcc


    ! this is be8(a,g)c12
    dd    = 130.0d0 * t9i32 * exp(-3.3364*t9i)
    ddd   = dd*(-1.5d0*t9i + 3.3364*t9i2)

    ee    = 2.510d+07 * t9i23 * exp(-23.57*t9i13 - t92*q2)
    dee   = ee*(-twoth*t9i + oneth*23.57*t9i43 - 2.0d0*t9*q2)
    
    ff    = 1.0d0 + 0.018*t913 + 5.249*t923 + 0.650*t9 + &
         19.176*t943 + 6.034*t953
    dff   = oneth*0.018*t9i23 + twoth*5.249*t9i13 + 0.650 &
         + fourth*19.176*t913 + fiveth*6.034*t923

    rbeac    = dd + ee * ff
    drbeacdt = ddd + dee * ff + ee * dff


    ! a factor
    xx    = rc28 * 1.35d-07 * t9i32 * exp(-24.811*t9i)
    dxx   = xx*(-1.5d0*t9i + 24.811*t9i2)


    ! high temperature rate
    if (t9.gt.0.08) then
       term    = 2.90d-16 * r2abe * rbeac + xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac &
                 + 2.90d-16 * r2abe * drbeacdt &
                 + dxx


    ! low temperature rate
    else
       uu   = 0.8d0*exp(-(0.025*t9i)**3.263)
       yy   = 0.2d0 + uu
       ! yy   = 0.01 + 0.2d0 + uu
       dyy  = uu * 3.263*(0.025*t9i)**2.263 * (0.025*t9i2)
       vv   = 4.0d0*exp(-(t9/0.025)**9.227)
       zz   = 1.0d0 + vv
       dzz  = vv * 9.227*(t9/0.025)**8.227 * 40.0d0
       aa   = 1.0d0/zz
       f1   = 0.01d0 + yy * aa
       ! f1   = yy * aa
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
    dfrdd = 2.0d0 * term * den

    rev    = 2.00d+20*t93*exp(-84.424*t9i)
    drevdt = rev*(3.0d0*t9i + 84.424*t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    drrdd = 0.0d0

    return
  end subroutine rate_tripalf


  subroutine rate_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                     aa,zz


    ! c12 + c12 reaction
    aa      = 1.0d0 + 0.0396*t9
    zz      = 1.0d0/aa

    t9a     = t9*zz
    dt9a    = (1.0d0 -  t9a*0.0396)*zz

    zz      = dt9a/t9a
    t9a13   = t9a**oneth
    dt9a13  = oneth*t9a13*zz
    
    t9a56   = t9a**fivsix
    dt9a56  = fivsix*t9a56*zz

    term    = 4.27d+26 * t9a56 * t9i32 * &
         exp(-84.165/t9a13 - 2.12d-03*t93)
    dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
         + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)

    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rr    = 0.0d0
    drrdt = 0.0d0
    drrdd = 0.0d0
    
    return
  end subroutine rate_c12c12


  subroutine rate_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23, &
                     t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz


    ! c12 + o16 reaction; see cf88 references 47-4
    if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
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
       term    = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                       + 106.594/t9a23*dt9a13 - zz*dcc)

    else
       ! term    = 2.6288035d-29
       term    = 0.0d0
       dtermdt = 0.0d0
    endif


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term

    rr    = 0.0d0
    drrdt = 0.0d0
    drrdd = 0.0d0
    
    return
  end subroutine rate_c12o16


  subroutine rate_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt


    ! o16 + o16
    term  = 7.10d36 * t9i23 * &
         exp(-135.93 * t9i13 - 0.629*t923 &
         - 0.445*t943 + 0.0103*t9*t9)

    dtermdt = -twoth*term*t9i &
         + term * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
         - fourth*0.445*t913 + 0.0206*t9)


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term

    rr    = 0.0d0
    drrdt = 0.0d0
    drrdd = 0.0d0
    
    return
  end subroutine rate_o16o16


  subroutine rate_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                     cc,dcc,term2,dterm2,rev,drevdt,q1
    parameter        (q1 = 1.0d0/2.515396d0)


    ! o16(a,g)ne20
    term1   = 9.37d9 * t9i23 * exp(-39.757*t9i13 - t92*q1)
    dterm1  = term1*(-twoth*t9i + oneth*39.757*t9i43 - 2.0d0*t9*q1)

    aa      = 62.1 * t9i32 * exp(-10.297*t9i)
    daa     = aa*(-1.5d0*t9i + 10.297*t9i2)
    
    bb      = 538.0d0 * t9i32 * exp(-12.226*t9i)
    dbb     = bb*(-1.5d0*t9i + 12.226*t9i2)

    cc      = 13.0d0 * t92 * exp(-20.093*t9i)
    dcc     = cc*(2.0d0*t9i + 20.093*t9i2)

    term2   = aa + bb + cc
    dterm2  = daa + dbb + dcc

    term    = term1 + term2
    dtermdt = dterm1 + dterm2


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term

    rev      = 5.65d+10*t932*exp(-54.937*t9i)
    drevdt   = rev*(1.5d0*t9i + 54.937*t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0

    return
  end subroutine rate_o16ag


  subroutine rate_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                     term2,dterm2,term3,dterm3,rev,drevdt,zz,rc102,q1
    parameter        (rc102 = 0.1d0, &
                      q1    = 1.0d0/4.923961d0)


    ! ne20(a,g)mg24
    aa   = 4.11d+11 * t9i23 * exp(-46.766*t9i13 - t92*q1)
    daa  = aa*(-twoth*t9i + oneth*46.766*t9i43 - 2.0d0*t9*q1)

    bb   = 1.0d0 + 0.009*t913 + 0.882*t923 + 0.055*t9 &
         + 0.749*t943 + 0.119*t953
    dbb  = oneth*0.009*t9i23 + twoth*0.882*t9i13 + 0.055 &
         + fourth*0.749*t913 + fiveth*0.119*t923
    
    term1  = aa * bb
    dterm1 = daa * bb + aa * dbb

    
    aa   = 5.27d+03 * t9i32 * exp(-15.869*t9i)
    daa  = aa*(-1.5d0*t9i + 15.869*t9i2)
    
    bb   = 6.51d+03 * t912 * exp(-16.223*t9i)
    dbb  = bb*(0.5d0*t9i + 16.223*t9i2)

    term2  = aa + bb
    dterm2 = daa + dbb


    aa   = 42.1 * t9i32 * exp(-9.115*t9i)
    daa  = aa*(-1.5d0*t9i + 9.115*t9i2)
    
    bb   =  32.0 * t9i23 * exp(-9.383*t9i)
    dbb  = bb*(-twoth*t9i + 9.383*t9i2)

    term3  = rc102 * (aa + bb)
    dterm3 = rc102 * (daa + dbb)


    aa  = 5.0d0*exp(-18.960*t9i)
    daa = aa*18.960*t9i2

    bb  = 1.0d0 + aa
    dbb = daa

    zz      = 1.0d0/bb
    term    = (term1 + term2 + term3)*zz
    dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 6.01d+10 * t932 * exp(-108.059*t9i)
    drevdt   = rev*(1.5d0*t9i + 108.059*t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0

    return
  end subroutine rate_ne20ag


  subroutine rate_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                     ff,dff,gg,dgg,hh,hhi,rev,drevdt,rc121
    parameter        (rc121 = 0.1d0)


    ! 24mg(a,g)28si
    aa    = 4.78d+01 * t9i32 * exp(-13.506*t9i)
    daa   = aa*(-1.5d0*t9i + 13.506*t9i2)

    bb    =  2.38d+03 * t9i32 * exp(-15.218*t9i)
    dbb   = bb*(-1.5d0*t9i + 15.218*t9i2)

    cc    = 2.47d+02 * t932 * exp(-15.147*t9i)
    dcc   = cc*(1.5d0*t9i + 15.147*t9i2)

    dd    = rc121 * 1.72d-09 * t9i32 * exp(-5.028*t9i)
    ddd   = dd*(-1.5d0*t9i + 5.028*t9i2)

    ee    = rc121* 1.25d-03 * t9i32 * exp(-7.929*t9i)
    dee   = ee*(-1.5d0*t9i + 7.929*t9i2)

    ff    = rc121 * 2.43d+01 * t9i * exp(-11.523*t9i)
    dff   = ff*(-t9i + 11.523*t9i2)

    gg    = 5.0d0*exp(-15.882*t9i)
    dgg   = gg*15.882*t9i2

    hh    = 1.0d0 + gg
    hhi   = 1.0d0/hh

    term    = (aa + bb + cc + dd + ee + ff) * hhi
    dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 6.27d+10 * t932 * exp(-115.862*t9i)
    drevdt   = rev*(1.5d0*t9i + 115.862*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_mg24ag


  subroutine rate_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                     ff,dff,gg,dgg,term1,dterm1,term2,dterm2, &
                     rev,drevdt,rc148,q1
    parameter        (rc148 = 0.1d0, &
                      q1    = 1.0d0/0.024649d0)


    ! 24mg(a,p)al27
    aa     = 1.10d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
    daa    = -twoth*aa*t9i + aa*(23.261*t9i43 - 2.0d0*t9*q1)

    bb     =  1.0d0 + 0.018*t913 + 12.85*t923 + 1.61*t9 &
         + 89.87*t943 + 28.66*t953
    dbb    = oneth*0.018*t9i23 + twoth*12.85*t9i13 + 1.61 &
           + fourth*89.87*t913 + fiveth*28.66*t923
    
    term1  = aa * bb
    dterm1 = daa * bb + aa * dbb

    aa     = 129.0d0 * t9i32 * exp(-2.517*t9i)
    daa    = -1.5d0*aa*t9i + aa*2.517*t9i2
    
    bb     = 5660.0d0 * t972 * exp(-3.421*t9i)
    dbb    = 3.5d0*bb*t9i +  bb*3.421*t9i2
    
    cc     = rc148 * 3.89d-08 * t9i32 * exp(-0.853*t9i)
    dcc    = -1.5d0*cc*t9i + cc*0.853*t9i2
    
    dd     = rc148 * 8.18d-09 * t9i32 * exp(-1.001*t9i)
    ddd    = -1.5d0*dd*t9i + dd*1.001*t9i2
    
    term2  = aa + bb + cc + dd
    dterm2 = daa + dbb + dcc + ddd

    ee     = oneth*exp(-9.792*t9i)
    dee    = ee*9.792*t9i2

    ff     =  twoth * exp(-11.773*t9i)
    dff    = ff*11.773*t9i2

    gg     = 1.0d0 + ee + ff
    dgg    = dee + dff

    term    = (term1 + term2)/gg
    dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


    ! the rates
    rev      = 1.81 * exp(-18.572*t9i)
    drevdt   = rev*18.572*t9i2

    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
    dfrdd = rev * term
    
    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    drrdd = term
    
    return
  end subroutine rate_mg24ap


  subroutine rate_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                      dd,ddd,ee,dee,ff,dff,gg,dgg


    ! al27(p,g)si28
    ! champagne 1996

    aa  = 1.32d+09 * t9i23 * exp(-23.26*t9i13)
    daa = aa*(-twoth*t9i + oneth*23.26*t9i43)

    bb  = 3.22d-10 * t9i32 * exp(-0.836*t9i)*0.17
    dbb = bb*(-1.5d0*t9i + 0.836*t9i2)

    cc  = 1.74d+00 * t9i32 * exp(-2.269*t9i)
    dcc = cc*(-1.5d0*t9i + 2.269*t9i2)

    dd  = 9.92d+00 * t9i32 * exp(-2.492*t9i)
    ddd = dd*(-1.5d0*t9i + 2.492*t9i2)

    ee  = 4.29d+01 * t9i32 * exp(-3.273*t9i)
    dee = ee*(-1.5d0*t9i + 3.273*t9i2)

    ff  = 1.34d+02 * t9i32 * exp(-3.654*t9i)
    dff = ff*(-1.5d0*t9i + 3.654*t9i2)

    gg  = 1.77d+04 * (t9**0.53) * exp(-4.588*t9i)
    dgg = gg*(0.53*t9i + 4.588*t9i2)

    term    = aa + bb + cc + dd + ee + ff + gg
    dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


    ! rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 1.13d+11 * t932 * exp(-134.434*t9i)
    drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

    rr    = rev * term
    drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_al27pg


  subroutine rate_al27pg_old(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                     ff,dff,gg,dgg,hh,dhh,xx,dxx,yy,dyy,zz,dzz,pp, &
                     rev,drevdt,rc147,q1
    parameter        (rc147 = 0.1d0, &
                      q1    = 1.0d0/0.024025d0)

    ! 27al(p,g)si28  cf88
    aa  = 1.67d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
    daa = aa*(-twoth*t9i + oneth*23.261*t9i43 - 2.0d0*t9*q1)

    bb  = 1.0d0 + 0.018*t913 + 5.81*t923 + 0.728*t9 &
         + 27.31*t943 + 8.71*t953
    dbb = oneth*0.018*t9i23 + twoth*5.81*t9i13 + 0.728 &
         + fourth*27.31*t913 + fiveth*8.71*t923
    
    cc  = aa*bb
    dcc = daa*bb + aa*dbb

    dd  = 2.20d+00 * t9i32 * exp(-2.269*t9i)
    ddd = dd*(-1.5d0*t9i + 2.269*t9i2)

    ee  = 1.22d+01 * t9i32 * exp(-2.491*t9i)
    dee = ee*(-1.5d0*t9i + 2.491*t9i2)

    ff  =  1.50d+04 * t9 * exp(-4.112*t9i)
    dff = ff*(t9i + 4.112*t9i2)

    gg  = rc147 * 6.50d-10 * t9i32 * exp(-0.853*t9i)
    dgg = gg*(-1.5d0*t9i + 0.853*t9i2)

    hh  = rc147 * 1.63d-10 * t9i32 * exp(-1.001*t9i)
    dhh = hh*(-1.5d0*t9i + 1.001*t9i2)

    xx     = oneth*exp(-9.792*t9i)
    dxx    = xx*9.792*t9i2

    yy     =  twoth * exp(-11.773*t9i)
    dyy    = yy*11.773*t9i2

    zz     = 1.0d0 + xx + yy
    dzz    = dxx + dyy

    pp      = 1.0d0/zz
    term    = (cc + dd + ee + ff + gg + hh)*pp
    dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term

    rev      = 1.13d+11*t932*exp(-134.434*t9i)
    drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0

    return
  end subroutine rate_al27pg_old


  subroutine rate_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


    ! si28(a,g)s32
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
    if (z .eq. 10.0) then
       daa = 0
    else
       daa   = 6.340d-2 + 2.0d0*2.541d-3*t9 - 3.0d0*2.900d-4*t92
    end if
    
    term    = 4.82d+22 * t9i23 * exp(-61.015 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 61.015*t9i13*(oneth*t9i*aa - daa))
      
    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 6.461d+10 * t932 * exp(-80.643*t9i)
    drevdt   = rev*(1.5d0*t9i + 80.643*t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0

    return
  end subroutine rate_si28ag


  subroutine rate_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


    ! si28(a,p)p31
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 2.798d-3 + 2.0d0*2.763d-3*t9 - 3.0d0*2.341d-4*t92
    end if
    
    term    = 4.16d+13 * t9i23 * exp(-25.631 * t9i13 * aa)
    dtermdt = -twoth*term*t9i + term*25.631*t9i13*(oneth*t9i*aa - daa)


    ! the rates
    rev      = 0.5825d0 * exp(-22.224*t9i)
    drevdt   = rev*22.224*t9i2
    
    fr    = den * rev * term
    dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
    dfrdd = rev * term
    
    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    drrdd = term
    
    return
  end subroutine rate_si28ap


  subroutine rate_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


    ! p31(p,g)s32
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.928d-1 - 2.0d0*1.540d-2*t9 + 3.0d0*6.444d-4*t92
    end if
    
    term    = 1.08d+16 * t9i23 * exp(-27.042 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 27.042*t9i13*(oneth*t9i*aa - daa))
    
    
    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 3.764d+10 * t932 * exp(-102.865*t9i)
    drevdt   = rev*(1.5d0*t9i + 102.865*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_p31pg


  subroutine rate_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


    ! s32(a,g)ar36
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 4.913d-2 + 2.0d0*4.637d-3*t9 - 3.0d0*4.067d-4*t92
    end if
    
    term    = 1.16d+24 * t9i23 * exp(-66.690 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 66.690*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 6.616d+10 * t932 * exp(-77.080*t9i)
    drevdt   = rev*(1.5d0*t9i + 77.080*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_s32ag
  

  subroutine rate_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


    ! s32(a,p)cl35
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
    if (z .eq. 10) then
       daa = 0.0d0
    else
       daa   = 1.041d-1 - 2.0d0*1.368d-2*t9 + 3.0d0*6.969d-4*t92
    end if
    
    term    = 1.27d+16 * t9i23 * exp(-31.044 * t9i13 * aa)
    dtermdt = -twoth*term*t9i + term*31.044*t9i13*(oneth*t9i*aa - daa)
    

    ! the rates
    rev      = 1.144 * exp(-21.643*t9i)
    drevdt   = rev*21.643*t9i2
    
    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    dfrdd = rev * term
    
    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    drrdd = term
    
    return
  end subroutine rate_s32ap


  subroutine rate_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt


    ! cl35(p,g)ar36
    aa    = 1.0d0 + 1.761d-1*t9 - 1.322d-2*t92 + 5.245d-4*t93
    daa   = 1.761d-1 - 2.0d0*1.322d-2*t9 + 3.0d0*5.245d-4*t92
      

    term    =  4.48d+16 * t9i23 * exp(-29.483 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 29.483*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 7.568d+10*t932*exp(-98.722*t9i)
    drevdt   = rev*(1.5d0*t9i + 98.722*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_cl35pg


  subroutine rate_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! ar36(a,g)ca40
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.458d-1 - 2.0d0*1.069d-2*t9 + 3.0d0*3.790d-4*t92
    end if
    
    term    = 2.81d+30 * t9i23 * exp(-78.271 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 78.271*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 6.740d+10 * t932 * exp(-81.711*t9i)
    drevdt   = rev*(1.5d0*t9i + 81.711*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_ar36ag


  subroutine rate_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


    ! ar36(a,p)k39
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 4.826d-3 - 2.0d0*5.534d-3*t9 + 3.0d0*4.021d-4*t92
    end if
    
    term    = 2.76d+13 * t9i23 * exp(-34.922 * t9i13 * aa)
    dtermdt = -twoth*term*t9i + term*34.922*t9i13*(oneth*t9i*aa - daa)


    ! the rates
    rev      = 1.128*exp(-14.959*t9i)
    drevdt   = rev*14.959*t9i2
    
    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    dfrdd = rev * term
    
    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    drrdd = term
    
    return
  end subroutine rate_ar36ap


  subroutine rate_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


    ! k39(p,g)ca40
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
    if (z .eq. 10) then
       daa = 0.0d0
    else
       daa   = 1.622d-1 - 2.0d0*1.119d-2*t9 + 3.0d0*3.910d-4*t92
    end if
    
    term    = 4.09d+16 * t9i23 * exp(-31.727 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 31.727*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 7.600d+10 * t932 * exp(-96.657*t9i)
    drevdt   = rev*(1.5d0*t9i + 96.657*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_k39pg


  subroutine rate_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    
    
    ! ca40(a,g)ti44
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.650d-2 + 2.0d0*5.973d-3*t9 - 3.0d0*3.889d-4*t92
    end if
    
    term    = 4.66d+24 * t9i23 * exp(-76.435 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 76.435*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 6.843d+10 * t932 * exp(-59.510*t9i)
    drevdt   = rev*(1.5d0*t9i + 59.510*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_ca40ag


  subroutine rate_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! ca40(a,p)sc43
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = -1.206d-2 + 2.0d0*7.753d-3*t9 - 3.0d0*5.071d-4*t92
    end if
    
    term    = 4.54d+14 * t9i23 * exp(-32.177 * t9i13 * aa)
    dtermdt = -twoth*term*t9i + term*32.177*t9i13*(oneth*t9i*aa - daa)


    ! the rates
    rev      = 2.229 * exp(-40.966*t9i)
    drevdt   = rev*40.966*t9i2
    
    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    dfrdd = rev * term
    
    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    drrdd = term
    
    return
  end subroutine rate_ca40ap


  subroutine rate_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    
    
    ! sc43(p,g)ca40
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.023d-1 - 2.0d0*2.242d-3*t9 - 3.0d0*5.463d-5*t92
    end if
    
    term    = 3.85d+16 * t9i23 * exp(-33.234 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 33.234*t9i13*(oneth*t9i*aa - daa))
    

    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 1.525d+11 * t932 * exp(-100.475*t9i)
    drevdt   = rev*(1.5d0*t9i + 100.475*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_sc43pg


  subroutine rate_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! ti44(a,g)cr48
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.066d-1 - 2.0d0*1.102d-2*t9 + 3.0d0*5.324d-4*t92
    end if
    
    term    = 1.37d+26 * t9i23 * exp(-81.227 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 81.227*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 6.928d+10*t932*exp(-89.289*t9i)
    drevdt   = rev*(1.5d0*t9i + 89.289*t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_ti44ag


  subroutine rate_ti44ap(temp,den, &
                         fr,dfrdt,dfrdd, &
                         rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! ti44(a,p)v47
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 2.655d-2 - 2.0d0*3.947d-3*t9 + 3.0d0*2.522d-4*t92
    end if
    
    term    = 6.54d+20 * t9i23 * exp(-66.678 * t9i13 * aa)
    dtermdt = -twoth*term*t9i + term*66.678*t9i13*(oneth*t9i*aa - daa)
    
    
    ! the rates
    rev      = 1.104 * exp(-4.723*t9i)
    drevdt   = rev*4.723*t9i2
    
    fr    = den * rev * term
    dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    dfrdd = rev * term
    
    rr    = den * term
    drrdt = den * dtermdt * 1.0d-9
    drrdd = term
    
    return
  end subroutine rate_ti44ap


  subroutine rate_v47pg(temp,den, &
                        fr,dfrdt,dfrdd, &
                        rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! v47(p,g)cr48
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 9.979d-2 - 2.0d0*2.269d-3*t9 - 3.0d0*6.662d-5*t92
    end if
    
    term    = 2.05d+17 * t9i23 * exp(-35.568 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 35.568*t9i13*(oneth*t9i*aa - daa))
    
    
    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 7.649d+10*t932*exp(-93.999*t9i)
    drevdt   = rev*(1.5d0*t9i + 93.999*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_v47pg


  subroutine rate_cr48ag(temp,den, &
                         fr,dfrdt,dfrdd, &
                         rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    
    
    ! cr48(a,g)fe52
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 6.325d-2 - 2.0d0*5.671d-3*t9 + 3.0d0*2.848d-4*t92
    end if
    
    term    = 1.04d+23 * t9i23 * exp(-81.420 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 81.420*t9i13*(oneth*t9i*aa - daa))
    
    
    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 7.001d+10 * t932 * exp(-92.177*t9i)
    drevdt   = rev*(1.5d0*t9i + 92.177*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_cr48ag


  subroutine rate_cr48ap(temp,den, &
                         fr,dfrdt,dfrdd, &
                         rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! cr48(a,p)mn51
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.384d-2 + 2.0d0*1.081d-3*t9 - 3.0d0*5.933d-5*t92
    end if
    
    term    = 1.83d+26 * t9i23 * exp(-86.741 * t9i13 * aa)
    dtermdt = -twoth*term*t9i + term*86.741*t9i13*(oneth*t9i*aa - daa)
    

    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 0.6087*exp(-6.510*t9i)
    drevdt   = rev*6.510*t9i2
    
    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    drrdd = rev * term
    
    return
  end subroutine rate_cr48ap


  subroutine rate_mn51pg(temp,den, &
                         fr,dfrdt,dfrdd, &
                         rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! mn51(p,g)fe52
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 8.922d-2 - 2.0d0*1.256d-3*t9 - 3.0d0*9.453d-5*t92
    end if
    
    term    = 3.77d+17 * t9i23 * exp(-37.516 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 37.516*t9i13*(oneth*t9i*aa - daa))
    

    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 1.150d+11*t932*exp(-85.667*t9i)
    drevdt   = rev*(1.5d0*t9i + 85.667*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_mn51pg


  subroutine rate_fe52ag(temp,den, &
                         fr,dfrdt,dfrdd, &
                         rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    

    ! fe52(a,g)ni56
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 7.846d-2 - 2.0d0*7.430d-3*t9 + 3.0d0*3.723d-4*t92
    end if
    
    term    = 1.05d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 91.674*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 7.064d+10*t932*exp(-92.850*t9i)
    drevdt   = rev*(1.5d0*t9i + 92.850*t9i2)

    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_fe52ag


  subroutine rate_fe52ap(temp,den, &
                         fr,dfrdt,dfrdd, &
                         rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd

    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    
    
    ! fe52(a,p)co55
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 1.367d-2 + 2.0d0*7.428d-4*t9 - 3.0d0*3.050d-5*t92
    end if
    
    term    = 1.30d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
    dtermdt = -twoth*term*t9i + term*91.674*t9i13*(oneth*t9i*aa - daa)


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 0.4597*exp(-9.470*t9i)
    drevdt   = rev*9.470*t9i2
    
    rr    = den * rev * term
    drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
    drrdd = rev * term
    
    return
  end subroutine rate_fe52ap
  

  subroutine rate_co55pg(temp,den, &
                         fr,dfrdt,dfrdd, &
                         rr,drrdt,drrdd)

    ! declare the pass
    double precision temp,den, &
                     fr,dfrdt,dfrdd, &
                     rr,drrdt,drrdd
    
    ! locals
    double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3
    
    
    ! co55(p,g)ni56
    z     = min(t9,10.0d0)
    z2    = z*z
    z3    = z2*z
    aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
    if (z .eq. 10.0) then
       daa = 0.0d0
    else
       daa   = 9.894d-2 - 2.0d0*3.131d-3*t9 - 3.0d0*2.160d-5*t92
    end if
    
    term    = 1.21d+18 * t9i23 * exp(-39.604 * t9i13 * aa)
    dtermdt = term*(-twoth*t9i + 39.604*t9i13*(oneth*t9i*aa - daa))


    ! the rates
    fr    = den * term
    dfrdt = den * dtermdt * 1.0d-9
    dfrdd = term
    
    rev      = 1.537d+11*t932*exp(-83.382*t9i)
    drevdt   = rev*(1.5d0*t9i + 83.382*t9i2)
    
    rr    = rev * term
    drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
    drrdd = 0.0d0
    
    return
  end subroutine rate_co55pg

