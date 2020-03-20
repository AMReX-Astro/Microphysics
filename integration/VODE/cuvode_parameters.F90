module cuvode_parameters_module

  use burn_type_module, only: neqs
#ifdef SIMPLIFIED_SDC
  use sdc_type_module, only : SVAR_EVOLVE
#endif
  use network, only : nspec

  implicit none

#ifdef TRUE_SDC
  integer, parameter :: VODE_NEQS = nspec + 2
#elif SIMPLIFIED_SDC
  integer, parameter :: VODE_NEQS = SVAR_EVOLVE
#else
  integer, parameter :: VODE_NEQS = neqs
#endif

end module cuvode_parameters_module
