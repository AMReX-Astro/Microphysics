! Various temperature factors -- these were originally defined in a common
! block in tfactors.dek.  Here we create a derived type to contain them 
! all -- they are typically passed into the rate routines.  
!
! Note: we comment out the ones that we don't explicitly use to save on
! computation.

module tfactors_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)          oneth,twoth,fourth,fiveth,elvnth,fivfour,onesix, &
                    fivsix,sevsix,onefif,sixfif,onesev,twosev, &
                    foursev

  parameter        (oneth   = 1.0e0_rt/3.0e0_rt, &
                    twoth   = 2.0e0_rt/3.0e0_rt, &
                    fourth  = 4.0e0_rt/3.0e0_rt, &
                    fiveth  = 5.0e0_rt/3.0e0_rt, &
                    elvnth  = 11.0e0_rt/3.0e0_rt, &
                    fivfour = 1.25e0_rt, &
                    onesix  = 1.0e0_rt/6.0e0_rt, &
                    fivsix  = 5.0e0_rt/6.0e0_rt, &
                    sevsix  = 7.0e0_rt/6.0e0_rt, &
                    onefif  = 0.2e0_rt, &
                    sixfif  = 1.2e0_rt, &
                    onesev  = 1.0e0_rt/7.0e0_rt, &
                    twosev  = 2.0e0_rt/7.0e0_rt, &
                    foursev = 4.0e0_rt/7.0e0_rt)


  type tf_t
     real(rt)         :: temp
     real(rt)         :: t9
     real(rt)         :: t92
     real(rt)         :: t93
     !real(rt)         :: t94
     real(rt)         :: t95
     !real(rt)         :: t96
     real(rt)         :: t912
     real(rt)         :: t932
     real(rt)         :: t952
     real(rt)         :: t972
     real(rt)         :: t913
     real(rt)         :: t923
     real(rt)         :: t943
     real(rt)         :: t953
     !real(rt)         :: t973
     !real(rt)         :: t9113
     !real(rt)         :: t914
     !real(rt)         :: t934
     !real(rt)         :: t954
     !real(rt)         :: t974
     !real(rt)         :: t915
     !real(rt)         :: t935
     !real(rt)         :: t945
     !real(rt)         :: t965
     real(rt)         :: t917
     !real(rt)         :: t927
     !real(rt)         :: t947
     !real(rt)         :: t918
     !real(rt)         :: t938
     !real(rt)         :: t958
     real(rt)         :: t9i
     real(rt)         :: t9i2
     real(rt)         :: t9i3
     real(rt)         :: t9i12
     real(rt)         :: t9i32
     !real(rt)         :: t9i52
     !real(rt)         :: t9i72
     real(rt)         :: t9i13 
     real(rt)         :: t9i23
     real(rt)         :: t9i43
     real(rt)         :: t9i53
     !real(rt)         :: t9i14
     !real(rt)         :: t9i34
     !real(rt)         :: t9i54
     !real(rt)         :: t9i15
     !real(rt)         :: t9i35 
     !real(rt)         :: t9i45
     !real(rt)         :: t9i65
     real(rt)         :: t9i17
     real(rt)         :: t9i27
     !real(rt)         :: t9i47
     !real(rt)         :: t9i18
     !real(rt)         :: t9i38
     !real(rt)         :: t9i58 
     !real(rt)         :: t916
     !real(rt)         :: t976
     !real(rt)         :: t9i76

  end type tf_t

contains

  subroutine get_tfactors(temp, tf)

    implicit none
    
    ! sets various popular temperature factors into common block this
    ! routine must be called before any of the rates are called

    ! declare the pass
    real(rt)         temp

    type(tf_t) :: tf

    !$gpu

    tf%temp = temp

    tf%t9    = temp * 1.0e-9_rt
    tf%t92   = tf%t9*tf%t9
    tf%t93   = tf%t9*tf%t92
    !tf%t94   = tf%t9*tf%t93
    !tf%t95   = tf%t9*tf%t94
    tf%t95   = tf%t92*tf%t93
    !tf%t96   = tf%t9*tf%t95
    
    tf%t912  = sqrt(tf%t9)
    tf%t932  = tf%t9*tf%t912
    tf%t952  = tf%t9*tf%t932
    !tf%t972  = tf%t9*tf%t952
    tf%t972  = tf%t92*tf%t932
    
    tf%t913  = tf%t9**oneth
    tf%t923  = tf%t913*tf%t913
    tf%t943  = tf%t9*tf%t913
    tf%t953  = tf%t9*tf%t923
    !tf%t973  = tf%t953*tf%t923
    !tf%t9113 = tf%t973*tf%t943

    !tf%t914  = tf%t9**(0.25e0_rt)
    !tf%t934  = tf%t914*tf%t914*tf%t914
    !tf%t954  = tf%t9*tf%t914
    !tf%t974  = tf%t9*tf%t934
    
    !tf%t915  = tf%t9**onefif
    !tf%t935  = tf%t915*tf%t915*tf%t915
    !tf%t945  = tf%t915 * tf%t935
    !tf%t965  = tf%t9 * tf%t915
    
    !tf%t916  = tf%t9**onesix
    !tf%t976  = tf%t9 * tf%t916
    !tf%t9i76 = 1.0e0_rt/tf%t976
    
    tf%t917  = tf%t9**onesev
    !tf%t927  = tf%t917*tf%t917
    !tf%t947  = tf%t927*tf%t927
    
    !tf%t918  = sqrt(tf%t914)
    !tf%t938  = tf%t918*tf%t918*tf%t918
    !tf%t958  = tf%t938*tf%t918*tf%t918
    
    tf%t9i   = 1.0e0_rt/tf%t9
    tf%t9i2  = tf%t9i*tf%t9i
    tf%t9i3  = tf%t9i2*tf%t9i
    
    tf%t9i12 = 1.0e0_rt/tf%t912
    tf%t9i32 = tf%t9i*tf%t9i12
    !tf%t9i52 = tf%t9i*tf%t9i32
    !tf%t9i72 = tf%t9i*tf%t9i52
    
    tf%t9i13 = 1.0e0_rt/tf%t913
    tf%t9i23 = tf%t9i13*tf%t9i13
    tf%t9i43 = tf%t9i*tf%t9i13
    tf%t9i53 = tf%t9i*tf%t9i23
    
    !tf%t9i14 = 1.0e0_rt/tf%t914
    !tf%t9i34 = tf%t9i14*tf%t9i14*tf%t9i14
    !tf%t9i54 = tf%t9i*tf%t9i14
    
    !tf%t9i15 = 1.0e0_rt/tf%t915
    !tf%t9i35 = tf%t9i15*tf%t9i15*tf%t9i15
    !tf%t9i45 = tf%t9i15 * tf%t9i35
    !tf%t9i65 = tf%t9i*tf%t9i15
    
    tf%t9i17 = 1.0e0_rt/tf%t917
    tf%t9i27 = tf%t9i17*tf%t9i17
    !tf%t9i47 = tf%t9i27*tf%t9i27
    
    !tf%t9i18 = 1.0e0_rt/tf%t918
    !tf%t9i38 = tf%t9i18*tf%t9i18*tf%t9i18
    !tf%t9i58 = tf%t9i38*tf%t9i18*tf%t9i18
    
  end subroutine get_tfactors

end module tfactors_module
