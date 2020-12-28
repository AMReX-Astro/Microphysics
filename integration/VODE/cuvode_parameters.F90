module cuvode_parameters_module

  use burn_type_module, only: neqs
  use network, only : nspec

  implicit none

#ifdef TRUE_SDC
  integer, parameter :: VODE_NEQS = nspec + 2
#else
  integer, parameter :: VODE_NEQS = neqs
#endif

end module cuvode_parameters_module
