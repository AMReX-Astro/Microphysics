module rate_type_module

  use actual_network, only : nrates

  implicit none

  integer, parameter :: num_rate_groups = 4

  type :: rate_t
     double precision :: rates(num_rate_groups, nrates)
  end type rate_t

end module rate_type_module
