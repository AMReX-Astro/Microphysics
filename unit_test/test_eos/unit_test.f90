subroutine init_unit_test(name, namlen) bind(C, name="init_unit_test")

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  call runtime_init(name, namlen)

end subroutine init_unit_test
