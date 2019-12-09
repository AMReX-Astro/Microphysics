module tfactors_module

  use microphysics_type_module, only: rt, THIRD

  implicit none

  type temp_t
     real(rt) :: t9
     real(rt) :: t9i
     real(rt) :: t9i13
     real(rt) :: t913
     real(rt) :: t953
     real(rt) :: lnt9
  end type temp_t

contains

  subroutine calc_tfactors(t9, tfactors)

    real(rt), intent(in   ) :: t9
    type (temp_t), intent(out) :: tfactors

    !$gpu

    tfactors%t9 = t9
    tfactors%t9i = 1.e0_rt / tfactors%t9
    tfactors%t9i13 = tfactors%t9i**THIRD
    tfactors%t913 = tfactors%t9**THIRD
    tfactors%t953 = tfactors%t9 * tfactors%t913 * tfactors%t913
    tfactors%lnt9 = log(tfactors%t9)

  end subroutine calc_tfactors

end module tfactors_module
