! Various temperature factors -- these were originally defined in a common
! block in tfactors.dek.  Here we create a derived type to contain them 
! all -- they are typically passed into the rate routines.  
!
! Note: we comment out the ones that we don't explicitly use to save on
! computation.

module tfactors_module

  implicit none

  double precision  oneth,twoth,fourth,fiveth,elvnth,fivfour,onesix, &
                    fivsix,sevsix,onefif,sixfif,onesev,twosev, &
                    foursev

  parameter        (oneth   = 1.0d0/3.0d0, &
                    twoth   = 2.0d0/3.0d0, &
                    fourth  = 4.0d0/3.0d0, &
                    fiveth  = 5.0d0/3.0d0, &
                    elvnth  = 11.0d0/3.0d0, &
                    fivfour = 1.25d0, &
                    onesix  = 1.0d0/6.0d0, &
                    fivsix  = 5.0d0/6.0d0, &
                    sevsix  = 7.0d0/6.0d0, &
                    onefif  = 0.2d0, &
                    sixfif  = 1.2d0, &
                    onesev  = 1.0d0/7.0d0, &
                    twosev  = 2.0d0/7.0d0, &
                    foursev = 4.0d0/7.0d0)


  type tf_t
     double precision :: temp
     double precision :: t9
     double precision :: t92
     double precision :: t93
     !double precision :: t94
     double precision :: t95
     !double precision :: t96
     double precision :: t912
     double precision :: t932
     double precision :: t952
     double precision :: t972
     double precision :: t913
     double precision :: t923
     double precision :: t943
     double precision :: t953
     !double precision :: t973
     !double precision :: t9113
     !double precision :: t914
     !double precision :: t934
     !double precision :: t954
     !double precision :: t974
     !double precision :: t915
     !double precision :: t935
     !double precision :: t945
     !double precision :: t965
     double precision :: t917
     !double precision :: t927
     !double precision :: t947
     !double precision :: t918
     !double precision :: t938
     !double precision :: t958
     double precision :: t9i
     double precision :: t9i2
     double precision :: t9i3
     double precision :: t9i12
     double precision :: t9i32
     !double precision :: t9i52
     !double precision :: t9i72
     double precision :: t9i13 
     double precision :: t9i23
     double precision :: t9i43
     double precision :: t9i53
     !double precision :: t9i14
     !double precision :: t9i34
     !double precision :: t9i54
     !double precision :: t9i15
     !double precision :: t9i35 
     !double precision :: t9i45
     !double precision :: t9i65
     double precision :: t9i17
     double precision :: t9i27
     !double precision :: t9i47
     !double precision :: t9i18
     !double precision :: t9i38
     !double precision :: t9i58 
     !double precision :: t916
     !double precision :: t976
     !double precision :: t9i76

  end type tf_t

contains

  subroutine get_tfactors(temp, tf)

    implicit none
    
    ! sets various popular temperature factors into common block this
    ! routine must be called before any of the rates are called

    ! declare the pass
    double precision temp

    type(tf_t) :: tf

    !$gpu

    tf%temp = temp

    tf%t9    = temp * 1.0d-9
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

    !tf%t914  = tf%t9**(0.25d0)
    !tf%t934  = tf%t914*tf%t914*tf%t914
    !tf%t954  = tf%t9*tf%t914
    !tf%t974  = tf%t9*tf%t934
    
    !tf%t915  = tf%t9**onefif
    !tf%t935  = tf%t915*tf%t915*tf%t915
    !tf%t945  = tf%t915 * tf%t935
    !tf%t965  = tf%t9 * tf%t915
    
    !tf%t916  = tf%t9**onesix
    !tf%t976  = tf%t9 * tf%t916
    !tf%t9i76 = 1.0d0/tf%t976
    
    tf%t917  = tf%t9**onesev
    !tf%t927  = tf%t917*tf%t917
    !tf%t947  = tf%t927*tf%t927
    
    !tf%t918  = sqrt(tf%t914)
    !tf%t938  = tf%t918*tf%t918*tf%t918
    !tf%t958  = tf%t938*tf%t918*tf%t918
    
    tf%t9i   = 1.0d0/tf%t9
    tf%t9i2  = tf%t9i*tf%t9i
    tf%t9i3  = tf%t9i2*tf%t9i
    
    tf%t9i12 = 1.0d0/tf%t912
    tf%t9i32 = tf%t9i*tf%t9i12
    !tf%t9i52 = tf%t9i*tf%t9i32
    !tf%t9i72 = tf%t9i*tf%t9i52
    
    tf%t9i13 = 1.0d0/tf%t913
    tf%t9i23 = tf%t9i13*tf%t9i13
    tf%t9i43 = tf%t9i*tf%t9i13
    tf%t9i53 = tf%t9i*tf%t9i23
    
    !tf%t9i14 = 1.0d0/tf%t914
    !tf%t9i34 = tf%t9i14*tf%t9i14*tf%t9i14
    !tf%t9i54 = tf%t9i*tf%t9i14
    
    !tf%t9i15 = 1.0d0/tf%t915
    !tf%t9i35 = tf%t9i15*tf%t9i15*tf%t9i15
    !tf%t9i45 = tf%t9i15 * tf%t9i35
    !tf%t9i65 = tf%t9i*tf%t9i15
    
    tf%t9i17 = 1.0d0/tf%t917
    tf%t9i27 = tf%t9i17*tf%t9i17
    !tf%t9i47 = tf%t9i27*tf%t9i27
    
    !tf%t9i18 = 1.0d0/tf%t918
    !tf%t9i38 = tf%t9i18*tf%t9i18*tf%t9i18
    !tf%t9i58 = tf%t9i38*tf%t9i18*tf%t9i18
    
  end subroutine get_tfactors

end module tfactors_module
