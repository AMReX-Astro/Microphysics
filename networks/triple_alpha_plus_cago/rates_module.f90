! This module contains the make_rates routine which calculates the
! rates for the triple alpha reaction network using the analytic
! expressions from caughlan and fowler (cf88).  We calculate both the
! forward, 3He4 --> C12, and backward, C12 --> 3He4 rates.  In
! addition, we include the C12 + He4 --> O16 forward and reverse rates
! also from cf88
!
! The rates calculated here are NOT screened and therefore an explicit 
! screening routine must be applied via screen.f
!


module rates_module

  use bl_types
  use bl_constants_module
  use network

  implicit none

contains

  subroutine make_rates(temp, dens, rates, dratesdt)

    !$acc routine seq
    
    ! rates given in terms of molar fractions

    real(kind=dp_t), intent(IN   ) :: temp, dens
    real(kind=dp_t), intent(  OUT) :: rates(nrates),dratesdt(nrates)
    
    real(kind=dp_t) :: t9r, t9r32, t9ri, t9ri2, t9, t9i, t913, t9i13, t923, &
                       t9i23, t943, t9i43, t932, t9i32, t953, t9i53, t92, t9i2

    real(kind=dp_t) ::  term,  dtermdt
    real(kind=dp_t) :: r2abe, dr2abedt
    real(kind=dp_t) :: rbeac, drbeacdt

    ! some scratch terms
    real(kind=dp_t) :: a, b, b1, b2, c
    real(kind=dp_t) :: dadt, dbdt, db1dt, db2dt, dcdt

    real(kind=dp_t), PARAMETER :: &
                     ONE_TWELVTH  = ONE / TWELVE, &
                     FIVE_SIXTHS  = FIVE * SIXTH, &
                     FIVE_THIRDS  = FIVE * THIRD, &
                     THREE_HALVES = THREE * HALF, &
                     T2T9         = 1.0e-9_dp_t


    t9r   = temp * T2T9
    t9r32 = t9r**THREE_HALVES
    t9ri  = ONE / t9r
    t9ri2 = t9ri**TWO
    t9    = min(t9r, 10.0e0_dp_t)
    t9i   = ONE / t9
    t913  = t9**THIRD
    t9i13 = ONE / t913
    t923  = t9**TWO3RD
    t9i23 = ONE / t923
    t943  = T9**FOUR3RD
    t9i43 = ONE / t943
    t932  = t9**THREE_HALVES
    t9i32 = ONE / t932
    t953  = t9**FIVE_THIRDS
    t9i53 = ONE / t953
    t92   = t9**TWO
    t9i2  = ONE / t92

    rates(:)    = ZERO
    dratesdt(:) = ZERO

    ! triple alpha to c12 (divided by 3! = 6)
    ! from cf88

    ! q = -0.092;     2 He4 --> Be8
    a    = (7.4e5_dp_t * t9i32) * dexp(-1.0663_dp_t * t9i)
    dadt = -a * THREE_HALVES * t9i + a * t9i2 * 1.0663_dp_t

    b    =  4.164e9_dp_t * t9i23 * dexp(-13.49_dp_t * t9i13 -                 &
                                         t92 / 9.604e-3_dp_t)
    dbdt = -b * TWO3RD * t9i + b * (13.49_dp_t * THIRD * t9i43 -              &
                                    TWO * t9 / 9.604e-3_dp_t)

    c    = ONE + 3.1e-2_dp_t * t913 + 8.009_dp_t * t923 +                     &
           1.732_dp_t * t9 + 49.883_dp_t * t943 + 27.426_dp_t * t953
    dcdt = 3.1e-2_dp_t * THIRD * t9i23 + 8.009_dp_t * TWO3RD * t9i13 +        &
           1.732_dp_t + 49.883_dp_t * FOUR3RD * t913 +                        &
           27.426_dp_t * FIVE_THIRDS * t923

    r2abe    = a + b * c
    dr2abedt = dadt + dbdt * c + b * dcdt


    ! q = 7.367;      He4 + Be8 --> C12
    a    = (130_dp_t * t9i32) * dexp(-3.3364_dp_t * t9i) 
    dadt = -a * THREE_HALVES * t9i + a * 3.3364_dp_t * t9i2
    
    b    = 2.51e7_dp_t * t9i23 * dexp(-23.57_dp_t * t9i13 -                   &
                                      t92 / 0.055225_dp_t)
    dbdt = b * TWO3RD * t9i + b * (23.57_dp_t * THIRD * t9i43 -               &
                                   TWO * t9 / 0.055225_dp_t)

    c    = ONE + 0.018_dp_t * t913 + 5.249_dp_t * t923 +                      &
           0.65_dp_t * t9 + 19.176_dp_t * t943 + 6.034_dp_t * t953
    dcdt = 0.018_dp_t * THIRD * t9i23 + 5.249_dp_t * TWO3RD * t9i13 +         &
           0.65_dp_t + 19.176_dp_t * FOUR3RD * t913 +                         &
           6.034_dp_t * FIVE_THIRDS * t923

    rbeac    = a + b * c
    drbeacdt = dadt + dbdt * c + b * dcdt

    ! q = 7.275;      total reaction

    a    = 2.9e-16_dp_t * r2abe * rbeac
    dadt = 2.9e-16_dp_t * (dr2abedt * rbeac + r2abe * drbeacdt)

    if (t9 .gt. 8e-2_dp_t) then
       b    = 0.1_dp_t * 1.35e-7_dp_t * t9i32 * dexp(-24.811_dp_t * t9i)
       dbdt = -b * TWO3RD * t9i + 24.811_dp_t * b * t9i2

       term    = a + b
       dtermdt = dadt + dbdt
    else
       b1    = ONE + FOUR * dexp(-(2.5e-2_dp_t * t9i)**3.263_dp_t)
       db1dt = FOUR * 3.263_dp_t * (2.5e-2_dp_t * t9i)**3.263_dp_t *          &
               t9i * dexp(-(2.5e-2_dp_t * t9i)**3.263_dp_t)

       b2    = ONE + FOUR * dexp(-(t9 / 2.5e-2_dp_t)**9.227_dp_t)
       db2dt = -FOUR * 9.227_dp_t * (t9 / 2.5e-2_dp_t)**9.227_dp_t *          &
               t9i * dexp(-(t9 / 2.5e-2_dp_t)**9.227_dp_t)

       b    = 1.e-2_dp_t + 0.2_dp_t * b1 / b2
       dbdt = 0.2_dp_t * (db1dt / b2 - b1 * db2dt / (b2 * b2))
            
       c    = 0.1_dp_t * 1.35e-7_dp_t * t9i32 * dexp(-24.811_dp_t * t9i)
       dcdt = -c * THREE_HALVES * t9i + 24.811_dp_t * c * t9i2

       term    = a * b + c
       dtermdt = dadt * b + a * dbdt + dcdt
    endif

    rates(ir3a)    = term * (dens * dens) / SIX
    dratesdt(ir3a) = dtermdt * T2T9 * (dens * dens) / SIX


    ! c12 + he4 to o16
    ! 1.7 time cf88 rate: see Weaver & Woosley Phy. Rep. 227 (1993)
    !                     and Garnett Nuc. Phys. A. 621  (1997)
    ! q = 7.162
    b1    = ONE + 0.0489_dp_t * t9i23
    db1dt = -0.0489_dp_t * TWO3RD * t9i53

    b2    = -32.120_dp_t * t9i13 - t92/(3.496_dp_t**2)
    db2dt = 32.120_dp_t * THIRD * t9i43 - TWO * t9 / (3.496_dp_t**2)

    a    = t9i2 * b1**2 * dexp(b2)
    dadt = a * (-TWO * t9i + TWO * db1dt / b1 + db2dt)
    
    !------------------------

    b2    = ONE + 0.2654_dp_t * t9i23
    db2dt = -0.2654_dp_t * TWO3RD * t9i53
    
    c    = -32.120_dp_t * t9i13
    dcdt = -THIRD * c * t9i

    b1    = t9i2 * b2**2 * dexp(c)
    db1dt = b1 * (-TWO * t9i + TWO * db2dt / b2 + dcdt)

    !------------------------

    c    = -27.499_dp_t * t9i
    dcdt = - c * t9i

    b2    = t9i32 * dexp(c)
    db2dt = b2 * (-THREE_HALVES * t9i + dcdt)

    !------------------------
    
    c    = t92 * t92 * t9 * dexp(-15.541_dp_t * t9i)
    dcdt = c * (FIVE * t9i + 15.541_dp_t * t9i2)

    !------------------------

    term    = 1.04e8_dp_t * a + 1.76e8_dp_t * b1 + &
              1.25e3_dp_t * b2 + 1.43e-2_dp_t * c
    dtermdt = 1.04e8_dp_t * dadt + 1.76e8_dp_t * db1dt + &
              1.25e3_dp_t * db2dt + 1.43e-2_dp_t * dcdt

    term    = 1.7_dp_t * term
    dtermdt = 1.7_dp_t * term

    rates(ircago)    = term * dens
    dratesdt(ircago) = dtermdt * T2T9 * dens

    return

  end subroutine make_rates

end module rates_module
