module rate_type_module

  use actual_network, only : nrates, num_rate_groups

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  type :: rate_t
     ! the temperature at which the rates were evaluated
     real(rt)         :: T_eval
     real(rt)         :: rates(num_rate_groups, nrates)
  end type rate_t

end module rate_type_module
