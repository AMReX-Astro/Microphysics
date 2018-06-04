
! This module is used to calculate reaction rates.
! Storage is used to handle common temperature factors used
! in the various reaction rates.
!
! This module is created via the make_rate_module.py routine.
module rates_module

  use amrex_constants_module

  implicit none

  type temp_t
     double precision :: t9
     double precision :: t9i
     double precision :: t9i13
     double precision :: t913
     double precision :: t953
     double precision :: lnt9
  end type temp_t

  double precision, parameter :: f17_to_o17 = exp(-4.53302d0)
  double precision, parameter :: ne18_to_f18 = exp(-0.880534d0)
  double precision, parameter :: ne19_to_f19 = exp(-3.21258d0)
  double precision, parameter :: o14_to_n14 = exp(-4.62412d0)
  double precision, parameter :: o15_to_n15 = exp(-5.1725d0)


contains

  function calc_tfactors(t9) result (tfactors)

    double precision, intent(in   ) :: t9
    type (temp_t) :: tfactors

    tfactors%t9 = t9
    tfactors%t9i = 1.d0 / tfactors%t9
    tfactors%t9i13 = tfactors%t9i**THIRD
    tfactors%t913 = tfactors%t9**THIRD
    tfactors%t953 = tfactors%t9 * tfactors%t913 * tfactors%t913
    tfactors%lnt9 = log(tfactors%t9)

  end function calc_tfactors


  subroutine rate_p_c12_to_n13(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1
    double precision :: dr0dt, dr1dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-3.77849d0)
    ct9i13 = tfactors%t9i13*(-5.10735d0)
    ct913 = tfactors%t913*(-2.24111d0)
    ct9 = tfactors%t9*(0.148883d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r0 = exp(17.5428d0 &
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
         +(-1.5d0))


    ct9i13 = tfactors%t9i13*(-13.692d0)
    ct913 = tfactors%t913*(-0.230881d0)
    ct9 = tfactors%t9*(4.44362d0)
    ct953 = tfactors%t953*(-3.15898d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r1 = exp(17.1482d0 &
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
         +(-0.666667d0))



    rate = r0 + r1 
    dratedt = dr0dt + dr1dt 

  end subroutine rate_p_c12_to_n13


  subroutine rate_f17_to_p_o16(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0
    double precision :: dr0dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-6.96583d0)
    ct9i13 = tfactors%t9i13*(-16.696d0)
    ct913 = tfactors%t913*(-1.16252d0)
    ct9 = tfactors%t9*(0.267703d0)
    ct953 = tfactors%t953*(-0.0338411d0)
    clnt9 = tfactors%lnt9*(0.833333d0)

    r0 = exp(40.9135d0 &
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
         +(0.833333d0))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_f17_to_p_o16


  subroutine rate_f17_to_o17(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt


    rate = f17_to_o17
    dratedt = ZERO

  end subroutine


  subroutine rate_p_f17_to_ne18(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1
    double precision :: dr0dt, dr1dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-0.0323504d0)
    ct9i13 = tfactors%t9i13*(-14.2191d0)
    ct913 = tfactors%t913*(34.0647d0)
    ct9 = tfactors%t9*(-16.5698d0)
    ct953 = tfactors%t953*(2.48116d0)
    clnt9 = tfactors%lnt9*(-2.13376d0)

    r0 = exp(-7.84708d0 &
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
         +(-2.13376d0))


    ct9i = tfactors%t9i*(-4.95969d0)
    ct9i13 = tfactors%t9i13*(-21.3249d0)
    ct913 = tfactors%t913*(-0.230774d0)
    ct9 = tfactors%t9*(0.917931d0)
    ct953 = tfactors%t953*(-0.0440377d0)
    clnt9 = tfactors%lnt9*(-7.36014d0)

    r1 = exp(27.5778d0 &
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
         +(-7.36014d0))



    rate = r0 + r1 
    dratedt = dr0dt + dr1dt 

  end subroutine rate_p_f17_to_ne18


  subroutine rate_he4_he4_he4_to_c12(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1, r2
    double precision :: dr0dt, dr1dt, dr2dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-1.02446d0)
    ct9i13 = tfactors%t9i13*(-23.57d0)
    ct913 = tfactors%t913*(20.4886d0)
    ct9 = tfactors%t9*(-12.9882d0)
    ct953 = tfactors%t953*(-20.0d0)
    clnt9 = tfactors%lnt9*(-2.16667d0)

    r0 = exp(-11.7884d0 &
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
         +(-2.16667d0))


    ct9i13 = tfactors%t9i13*(-37.06d0)
    ct913 = tfactors%t913*(29.3493d0)
    ct9 = tfactors%t9*(-115.507d0)
    ct953 = tfactors%t953*(-10.0d0)
    clnt9 = tfactors%lnt9*(-1.33333d0)

    r1 = exp(-0.971052d0 &
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
         +(-1.33333d0))


    ct9i = tfactors%t9i*(-4.12656d0)
    ct9i13 = tfactors%t9i13*(-13.49d0)
    ct913 = tfactors%t913*(21.4259d0)
    ct9 = tfactors%t9*(-1.34769d0)
    ct953 = tfactors%t953*(0.0879816d0)
    clnt9 = tfactors%lnt9*(-13.1653d0)

    r2 = exp(-24.3505d0 &
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
         +(-13.1653d0))



    rate = r0 + r1 + r2 
    dratedt = dr0dt + dr1dt + dr2dt 

  end subroutine rate_he4_he4_he4_to_c12


  subroutine rate_p_n14_to_o15(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1, r2, r3
    double precision :: dr0dt, dr1dt, dr2dt, dr3dt

    rate = ZERO
    dratedt = ZERO

    ct9i13 = tfactors%t9i13*(-15.193d0)
    ct913 = tfactors%t913*(-4.63975d0)
    ct9 = tfactors%t9*(9.73458d0)
    ct953 = tfactors%t953*(-9.55051d0)
    clnt9 = tfactors%lnt9*(0.333333d0)

    r0 = exp(20.1169d0 &
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
         +(0.333333d0))


    ct9i13 = tfactors%t9i13*(-15.193d0)
    ct913 = tfactors%t913*(-0.161954d0)
    ct9 = tfactors%t9*(-7.52123d0)
    ct953 = tfactors%t953*(-0.987565d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r1 = exp(17.01d0 &
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
         +(-0.666667d0))


    ct9i = tfactors%t9i*(-4.891d0)
    clnt9 = tfactors%lnt9*(0.0682d0)

    r2 = exp(6.73578d0 &
         +ct9i &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -ct9i &
         +(0.0682d0))


    ct9i = tfactors%t9i*(-2.998d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r3 = exp(7.65444d0 &
         +ct9i &
         +clnt9)
    dr3dt = r3 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))



    rate = r0 + r1 + r2 + r3 
    dratedt = dr0dt + dr1dt + dr2dt + dr3dt 

  end subroutine rate_p_n14_to_o15


  subroutine rate_he4_ne18_to_p_na21(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1, r2, r3, r4, r5, r6
    double precision :: dr0dt, dr1dt, dr2dt, dr3dt, dr4dt, dr5dt, dr6dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-24.7176d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r0 = exp(19.4058d0 &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-0.452576d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r1 = exp(-137.358d0 &
         +ct9i &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-10.885d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r2 = exp(6.39797d0 &
         +ct9i &
         +clnt9)
    dr2dt = r2 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-7.45009d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r3 = exp(-1.15641d0 &
         +ct9i &
         +clnt9)
    dr3dt = r3 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-5.97632d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r4 = exp(-6.65137d0 &
         +ct9i &
         +clnt9)
    dr4dt = r4 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-3.14078d0)
    ct9i13 = tfactors%t9i13*(-17.7759d0)
    ct913 = tfactors%t913*(36.0724d0)
    ct9 = tfactors%t9*(-5.34039d0)
    ct953 = tfactors%t953*(0.382679d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r5 = exp(-13.3392d0 &
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
         +(-1.5d0))


    ct9i = tfactors%t9i*(-2.81989d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r6 = exp(-28.6929d0 &
         +ct9i &
         +clnt9)
    dr6dt = r6 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))



    rate = r0 + r1 + r2 + r3 + r4 + r5 + r6 
    dratedt = dr0dt + dr1dt + dr2dt + dr3dt + dr4dt + dr5dt + dr6dt 

  end subroutine rate_he4_ne18_to_p_na21


  subroutine rate_ne18_to_f18(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt


    rate = ne18_to_f18
    dratedt = ZERO

  end subroutine


  subroutine rate_ne19_to_f19(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt


    rate = ne19_to_f19
    dratedt = ZERO

  end subroutine


  subroutine rate_p_ne19_to_na20(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1
    double precision :: dr0dt, dr1dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-5.07623d0)
    ct913 = tfactors%t913*(1.23704d0)
    ct9 = tfactors%t9*(0.337618d0)
    ct953 = tfactors%t953*(-0.0562825d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r0 = exp(5.63289d0 &
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
         +(-1.5d0))


    ct9i13 = tfactors%t9i13*(-19.5908d0)
    ct913 = tfactors%t913*(-2.37696d0)
    ct9 = tfactors%t9*(3.26815d0)
    ct953 = tfactors%t953*(-1.06524d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r1 = exp(17.822d0 &
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
         +(-0.666667d0))



    rate = r0 + r1 
    dratedt = dr0dt + dr1dt 

  end subroutine rate_p_ne19_to_na20


  subroutine rate_he4_o14_to_p_f17(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1, r2, r3, r4, r5
    double precision :: dr0dt, dr1dt, dr2dt, dr3dt, dr4dt, dr5dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-12.0223d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r0 = exp(12.1289d0 &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-26.0d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r1 = exp(18.6518d0 &
         +ct9i &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i13 = tfactors%t9i13*(-39.388d0)
    ct913 = tfactors%t913*(-17.4673d0)
    ct9 = tfactors%t9*(35.3029d0)
    ct953 = tfactors%t953*(-24.8162d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r2 = exp(40.8358d0 &
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
         +(-0.666667d0))


    ct9i = tfactors%t9i*(-22.51d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r3 = exp(16.3087d0 &
         +ct9i &
         +clnt9)
    dr3dt = r3 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-13.6d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r4 = exp(11.1184d0 &
         +ct9i &
         +clnt9)
    dr4dt = r4 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i = tfactors%t9i*(-0.453036d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r5 = exp(-106.091d0 &
         +ct9i &
         +clnt9)
    dr5dt = r5 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))



    rate = r0 + r1 + r2 + r3 + r4 + r5 
    dratedt = dr0dt + dr1dt + dr2dt + dr3dt + dr4dt + dr5dt 

  end subroutine rate_he4_o14_to_p_f17


  subroutine rate_o14_to_n14(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt


    rate = o14_to_n14
    dratedt = ZERO

  end subroutine


  subroutine rate_he4_o15_to_ne19(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1, r2
    double precision :: dr0dt, dr1dt, dr2dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-5.88439d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r0 = exp(-0.0452465d0 &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i13 = tfactors%t9i13*(-39.578d0)
    ct953 = tfactors%t953*(-3.0d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r1 = exp(26.2914d0 &
         +ct9i13 &
         +ct953 &
         +clnt9)
    dr1dt = r1 * tfactors%t9i * ( &
         -THIRD*ct9i13 &
         +FIVE3RD*ct953 &
         +(-0.666667d0))


    ct9i = tfactors%t9i*(-4.20439d0)
    ct9i13 = tfactors%t9i13*(-3.24609d0)
    ct913 = tfactors%t913*(44.4647d0)
    ct9 = tfactors%t9*(-9.79962d0)
    ct953 = tfactors%t953*(0.841782d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r2 = exp(-32.2496d0 &
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
         +(-1.5d0))



    rate = r0 + r1 + r2 
    dratedt = dr0dt + dr1dt + dr2dt 

  end subroutine rate_he4_o15_to_ne19


  subroutine rate_o15_to_n15(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt


    rate = o15_to_n15
    dratedt = ZERO

  end subroutine


  subroutine rate_he4_o16_to_ne20(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0, r1, r2
    double precision :: dr0dt, dr1dt, dr2dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-10.3585d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r0 = exp(3.88571d0 &
         +ct9i &
         +clnt9)
    dr0dt = r0 * tfactors%t9i * ( &
         -ct9i &
         +(-1.5d0))


    ct9i13 = tfactors%t9i13*(-39.7262d0)
    ct913 = tfactors%t913*(-0.210799d0)
    ct9 = tfactors%t9*(0.442879d0)
    ct953 = tfactors%t953*(-0.0797753d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r1 = exp(23.903d0 &
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
         +(-0.666667d0))


    ct9i = tfactors%t9i*(-12.7643d0)
    ct913 = tfactors%t913*(-3.65925d0)
    ct9 = tfactors%t9*(0.714224d0)
    ct953 = tfactors%t953*(-0.00107508d0)
    clnt9 = tfactors%lnt9*(-1.5d0)

    r2 = exp(9.50848d0 &
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
         +(-1.5d0))



    rate = r0 + r1 + r2 
    dratedt = dr0dt + dr1dt + dr2dt 

  end subroutine rate_he4_o16_to_ne20


  subroutine rate_p_o16_to_f17(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0
    double precision :: dr0dt

    rate = ZERO
    dratedt = ZERO

    ct9i13 = tfactors%t9i13*(-16.696d0)
    ct913 = tfactors%t913*(-1.16252d0)
    ct9 = tfactors%t9*(0.267703d0)
    ct953 = tfactors%t953*(-0.0338411d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r0 = exp(19.0904d0 &
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
         +(-0.666667d0))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_p_o16_to_f17


  subroutine rate_he4_si26_to_p_p29(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0
    double precision :: dr0dt

    rate = ZERO
    dratedt = ZERO

    ct9i13 = tfactors%t9i13*(-59.3013d0)
    ct913 = tfactors%t913*(0.480742d0)
    ct9 = tfactors%t9*(-0.834505d0)
    ct953 = tfactors%t953*(0.0621841d0)
    clnt9 = tfactors%lnt9*(-0.666667d0)

    r0 = exp(48.8732d0 &
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
         +(-0.666667d0))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_he4_si26_to_p_p29


  subroutine rate_he4_ti44_to_p_v47(tfactors,rate,dratedt)

    type (temp_t),    intent(in   ) :: tfactors
    double precision, intent(  out) :: rate
    double precision, intent(  out) :: dratedt

    double precision :: ct9i, ct9i13, ct913, ct9, ct953, clnt9

    double precision :: r0
    double precision :: dr0dt

    rate = ZERO
    dratedt = ZERO

    ct9i = tfactors%t9i*(-9.07869d0)
    ct9i13 = tfactors%t9i13*(5.56533d0)
    ct913 = tfactors%t913*(18.4415d0)
    ct9 = tfactors%t9*(-4.10095d0)
    ct953 = tfactors%t953*(0.24244d0)
    clnt9 = tfactors%lnt9*(16.0516d0)

    r0 = exp(-34.2468d0 &
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
         +(16.0516d0))



    rate = r0 
    dratedt = dr0dt 

  end subroutine rate_he4_ti44_to_p_v47



end module rates_module
