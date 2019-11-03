module tfactors_module

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

contains

  subroutine calc_tfactors(t9, tfactors)

    double precision, intent(in   ) :: t9
    type (temp_t), intent(out) :: tfactors

    !$gpu

    tfactors%t9 = t9
    tfactors%t9i = 1.d0 / tfactors%t9
    tfactors%t9i13 = tfactors%t9i**THIRD
    tfactors%t913 = tfactors%t9**THIRD
    tfactors%t953 = tfactors%t9 * tfactors%t913 * tfactors%t913
    tfactors%lnt9 = log(tfactors%t9)

  end subroutine calc_tfactors

end module tfactors_module
