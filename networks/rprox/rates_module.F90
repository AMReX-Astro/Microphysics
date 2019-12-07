
! This module is used to calculate reaction rates.
! Storage is used to handle common temperature factors used
! in the various reaction rates.
!
! This module is created via the make_rate_module.py routine.
module rates_module

  use amrex_constants_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  real(rt)        , parameter :: f17_to_o17 = exp(-4.53302e0_rt)
  real(rt)        , parameter :: ne18_to_f18 = exp(-0.880534e0_rt)
  real(rt)        , parameter :: ne19_to_f19 = exp(-3.21258e0_rt)
  real(rt)        , parameter :: o14_to_n14 = exp(-4.62412e0_rt)
  real(rt)        , parameter :: o15_to_n15 = exp(-5.1725e0_rt)


contains

  subroutine rate_p_c12_to_n13(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1
    real(rt)         :: dr0dt, dr1dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-3.77849e0_rt)
    ct9i13 = tfactors%t9i13*(-5.10735e0_rt)
    ct913 = tfactors%t913*(-2.24111e0_rt)
    ct9 = tfactors%t9*(0.148883e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r0 = exp(17.5428e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +(-1.5e0_rt))


    ct9i13 = tfactors%t9i13*(-13.692e0_rt)
    ct913 = tfactors%t913*(-0.230881e0_rt)
    ct9 = tfactors%t9*(4.44362e0_rt)
    ct953 = tfactors%t953*(-3.15898e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r1 = exp(17.1482e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))



    rate = r0 + r1 
    dratedt = dr0dt + dr1dt 

  end subroutine rate_p_c12_to_n13


  subroutine rate_f17_to_p_o16(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0
    real(rt)         :: dr0dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-6.96583e0_rt)
    ct9i13 = tfactors%t9i13*(-16.696e0_rt)
    ct913 = tfactors%t913*(-1.16252e0_rt)
    ct9 = tfactors%t9*(0.267703e0_rt)
    ct953 = tfactors%t953*(-0.0338411e0_rt)
    clnt9 = tfactors%lnt9*(0.833333e0_rt)

    r0 = exp(40.9135e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(0.833333e0_rt))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_f17_to_p_o16


  subroutine rate_f17_to_o17(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    !$gpu

    rate = f17_to_o17
    dratedt = ZERO

  end subroutine rate_f17_to_o17


  subroutine rate_p_f17_to_ne18(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1
    real(rt)         :: dr0dt, dr1dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-0.0323504e0_rt)
    ct9i13 = tfactors%t9i13*(-14.2191e0_rt)
    ct913 = tfactors%t913*(34.0647e0_rt)
    ct9 = tfactors%t9*(-16.5698e0_rt)
    ct953 = tfactors%t953*(2.48116e0_rt)
    clnt9 = tfactors%lnt9*(-2.13376e0_rt)

    r0 = exp(-7.84708e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-2.13376e0_rt))


    ct9i = tfactors%t9i*(-4.95969e0_rt)
    ct9i13 = tfactors%t9i13*(-21.3249e0_rt)
    ct913 = tfactors%t913*(-0.230774e0_rt)
    ct9 = tfactors%t9*(0.917931e0_rt)
    ct953 = tfactors%t953*(-0.0440377e0_rt)
    clnt9 = tfactors%lnt9*(-7.36014e0_rt)

    r1 = exp(27.5778e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-7.36014e0_rt))



    rate = r0 + r1 
    dratedt = dr0dt + dr1dt 

  end subroutine rate_p_f17_to_ne18


  subroutine rate_he4_he4_he4_to_c12(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1, r2
    real(rt)         :: dr0dt, dr1dt, dr2dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-1.02446e0_rt)
    ct9i13 = tfactors%t9i13*(-23.57e0_rt)
    ct913 = tfactors%t913*(20.4886e0_rt)
    ct9 = tfactors%t9*(-12.9882e0_rt)
    ct953 = tfactors%t953*(-20.0e0_rt)
    clnt9 = tfactors%lnt9*(-2.16667e0_rt)

    r0 = exp(-11.7884e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-2.16667e0_rt))


    ct9i13 = tfactors%t9i13*(-37.06e0_rt)
    ct913 = tfactors%t913*(29.3493e0_rt)
    ct9 = tfactors%t9*(-115.507e0_rt)
    ct953 = tfactors%t953*(-10.0e0_rt)
    clnt9 = tfactors%lnt9*(-1.33333e0_rt)

    r1 = exp(-0.971052e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-1.33333e0_rt))


    ct9i = tfactors%t9i*(-4.12656e0_rt)
    ct9i13 = tfactors%t9i13*(-13.49e0_rt)
    ct913 = tfactors%t913*(21.4259e0_rt)
    ct9 = tfactors%t9*(-1.34769e0_rt)
    ct953 = tfactors%t953*(0.0879816e0_rt)
    clnt9 = tfactors%lnt9*(-13.1653e0_rt)

    r2 = exp(-24.3505e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-13.1653e0_rt))



    rate = r0 + r1 + r2 
    dratedt = dr0dt + dr1dt + dr2dt 

  end subroutine rate_he4_he4_he4_to_c12


  subroutine rate_p_n14_to_o15(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1, r2, r3
    real(rt)         :: dr0dt, dr1dt, dr2dt, dr3dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i13 = tfactors%t9i13*(-15.193e0_rt)
    ct913 = tfactors%t913*(-4.63975e0_rt)
    ct9 = tfactors%t9*(9.73458e0_rt)
    ct953 = tfactors%t953*(-9.55051e0_rt)
    clnt9 = tfactors%lnt9*(0.333333e0_rt)

    r0 = exp(20.1169e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(0.333333e0_rt))


    ct9i13 = tfactors%t9i13*(-15.193e0_rt)
    ct913 = tfactors%t913*(-0.161954e0_rt)
    ct9 = tfactors%t9*(-7.52123e0_rt)
    ct953 = tfactors%t953*(-0.987565e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r1 = exp(17.01e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))


    ct9i = tfactors%t9i*(-4.891e0_rt)
    clnt9 = tfactors%lnt9*(0.0682e0_rt)

    r2 = exp(6.73578e0_rt &
         +ct9i &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -ct9i &
         +(0.0682e0_rt))


    ct9i = tfactors%t9i*(-2.998e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r3 = exp(7.65444e0_rt &
         +ct9i &
         +clnt9)
    dr3dt = r3 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))



    rate = r0 + r1 + r2 + r3 
    dratedt = dr0dt + dr1dt + dr2dt + dr3dt 

  end subroutine rate_p_n14_to_o15


  subroutine rate_he4_ne18_to_p_na21(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1, r2, r3, r4, r5, r6
    real(rt)         :: dr0dt, dr1dt, dr2dt, dr3dt, dr4dt, dr5dt, dr6dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-24.7176e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r0 = exp(19.4058e0_rt &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-0.452576e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r1 = exp(-137.358e0_rt &
         +ct9i &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-10.885e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r2 = exp(6.39797e0_rt &
         +ct9i &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-7.45009e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r3 = exp(-1.15641e0_rt &
         +ct9i &
         +clnt9)
    dr3dt = r3 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-5.97632e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r4 = exp(-6.65137e0_rt &
         +ct9i &
         +clnt9)
    dr4dt = r4 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-3.14078e0_rt)
    ct9i13 = tfactors%t9i13*(-17.7759e0_rt)
    ct913 = tfactors%t913*(36.0724e0_rt)
    ct9 = tfactors%t9*(-5.34039e0_rt)
    ct953 = tfactors%t953*(0.382679e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r5 = exp(-13.3392e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr5dt = r5 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-2.81989e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r6 = exp(-28.6929e0_rt &
         +ct9i &
         +clnt9)
    dr6dt = r6 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))



    rate = r0 + r1 + r2 + r3 + r4 + r5 + r6 
    dratedt = dr0dt + dr1dt + dr2dt + dr3dt + dr4dt + dr5dt + dr6dt 

  end subroutine rate_he4_ne18_to_p_na21


  subroutine rate_ne18_to_f18(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    !$gpu

    rate = ne18_to_f18
    dratedt = ZERO

  end subroutine rate_ne18_to_f18


  subroutine rate_ne19_to_f19(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    !$gpu

    rate = ne19_to_f19
    dratedt = ZERO

  end subroutine rate_ne19_to_f19


  subroutine rate_p_ne19_to_na20(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1
    real(rt)         :: dr0dt, dr1dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-5.07623e0_rt)
    ct913 = tfactors%t913*(1.23704e0_rt)
    ct9 = tfactors%t9*(0.337618e0_rt)
    ct953 = tfactors%t953*(-0.0562825e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r0 = exp(5.63289e0_rt &
         +ct9i &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-1.5e0_rt))


    ct9i13 = tfactors%t9i13*(-19.5908e0_rt)
    ct913 = tfactors%t913*(-2.37696e0_rt)
    ct9 = tfactors%t9*(3.26815e0_rt)
    ct953 = tfactors%t953*(-1.06524e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r1 = exp(17.822e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))



    rate = r0 + r1 
    dratedt = dr0dt + dr1dt 

  end subroutine rate_p_ne19_to_na20


  subroutine rate_he4_o14_to_p_f17(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1, r2, r3, r4, r5
    real(rt)         :: dr0dt, dr1dt, dr2dt, dr3dt, dr4dt, dr5dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-12.0223e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r0 = exp(12.1289e0_rt &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-26.0e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r1 = exp(18.6518e0_rt &
         +ct9i &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i13 = tfactors%t9i13*(-39.388e0_rt)
    ct913 = tfactors%t913*(-17.4673e0_rt)
    ct9 = tfactors%t9*(35.3029e0_rt)
    ct953 = tfactors%t953*(-24.8162e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r2 = exp(40.8358e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))


    ct9i = tfactors%t9i*(-22.51e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r3 = exp(16.3087e0_rt &
         +ct9i &
         +clnt9)
    dr3dt = r3 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-13.6e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r4 = exp(11.1184e0_rt &
         +ct9i &
         +clnt9)
    dr4dt = r4 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i = tfactors%t9i*(-0.453036e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r5 = exp(-106.091e0_rt &
         +ct9i &
         +clnt9)
    dr5dt = r5 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))



    rate = r0 + r1 + r2 + r3 + r4 + r5 
    dratedt = dr0dt + dr1dt + dr2dt + dr3dt + dr4dt + dr5dt 

  end subroutine rate_he4_o14_to_p_f17


  subroutine rate_o14_to_n14(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    !$gpu

    rate = o14_to_n14
    dratedt = ZERO

  end subroutine rate_o14_to_n14


  subroutine rate_he4_o15_to_ne19(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1, r2
    real(rt)         :: dr0dt, dr1dt, dr2dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-5.88439e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r0 = exp(-0.0452465e0_rt &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i13 = tfactors%t9i13*(-39.578e0_rt)
    ct953 = tfactors%t953*(-3.0e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r1 = exp(26.2914e0_rt &
         +ct9i13 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))


    ct9i = tfactors%t9i*(-4.20439e0_rt)
    ct9i13 = tfactors%t9i13*(-3.24609e0_rt)
    ct913 = tfactors%t913*(44.4647e0_rt)
    ct9 = tfactors%t9*(-9.79962e0_rt)
    ct953 = tfactors%t953*(0.841782e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r2 = exp(-32.2496e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-1.5e0_rt))



    rate = r0 + r1 + r2 
    dratedt = dr0dt + dr1dt + dr2dt 

  end subroutine rate_he4_o15_to_ne19


  subroutine rate_o15_to_n15(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    !$gpu

    rate = o15_to_n15
    dratedt = ZERO

  end subroutine rate_o15_to_n15


  subroutine rate_he4_o16_to_ne20(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0, r1, r2
    real(rt)         :: dr0dt, dr1dt, dr2dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-10.3585e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r0 = exp(3.88571e0_rt &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5e0_rt))


    ct9i13 = tfactors%t9i13*(-39.7262e0_rt)
    ct913 = tfactors%t913*(-0.210799e0_rt)
    ct9 = tfactors%t9*(0.442879e0_rt)
    ct953 = tfactors%t953*(-0.0797753e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r1 = exp(23.903e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))


    ct9i = tfactors%t9i*(-12.7643e0_rt)
    ct913 = tfactors%t913*(-3.65925e0_rt)
    ct9 = tfactors%t9*(0.714224e0_rt)
    ct953 = tfactors%t953*(-0.00107508e0_rt)
    clnt9 = tfactors%lnt9*(-1.5e0_rt)

    r2 = exp(9.50848e0_rt &
         +ct9i &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -ct9i &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-1.5e0_rt))



    rate = r0 + r1 + r2 
    dratedt = dr0dt + dr1dt + dr2dt 

  end subroutine rate_he4_o16_to_ne20


  subroutine rate_p_o16_to_f17(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0
    real(rt)         :: dr0dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i13 = tfactors%t9i13*(-16.696e0_rt)
    ct913 = tfactors%t913*(-1.16252e0_rt)
    ct9 = tfactors%t9*(0.267703e0_rt)
    ct953 = tfactors%t953*(-0.0338411e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r0 = exp(19.0904e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_p_o16_to_f17


  subroutine rate_he4_si26_to_p_p29(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0
    real(rt)         :: dr0dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i13 = tfactors%t9i13*(-59.3013e0_rt)
    ct913 = tfactors%t913*(0.480742e0_rt)
    ct9 = tfactors%t9*(-0.834505e0_rt)
    ct953 = tfactors%t953*(0.0621841e0_rt)
    clnt9 = tfactors%lnt9*(-0.666667e0_rt)

    r0 = exp(48.8732e0_rt &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(-0.666667e0_rt))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_he4_si26_to_p_p29


  subroutine rate_he4_ti44_to_p_v47(tfactors,rate,dratedt)

    use tfactors_module, only : temp_t

    use amrex_fort_module, only : rt => amrex_real
    type (temp_t),    intent(in   ) :: tfactors
    real(rt)        , intent(  out) :: rate
    real(rt)        , intent(  out) :: dratedt

    real(rt)         :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    real(rt)         :: r0
    real(rt)         :: dr0dt

    !$gpu

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-9.07869e0_rt)
    ct9i13 = tfactors%t9i13*(5.56533e0_rt)
    ct913 = tfactors%t913*(18.4415e0_rt)
    ct9 = tfactors%t9*(-4.10095e0_rt)
    ct953 = tfactors%t953*(0.24244e0_rt)
    clnt9 = tfactors%lnt9*(16.0516e0_rt)

    r0 = exp(-34.2468e0_rt &
         +ct9i &
         +ct9i13 &
         +ct913 &
         +ct9 &
         +ct953 &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         -THIRD*ct9i13 &
         +THIRD*ct913 &
         +ct9 &
         +FIVE3RD*ct953 &
         +(16.0516e0_rt))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_he4_ti44_to_p_v47

end module rates_module
