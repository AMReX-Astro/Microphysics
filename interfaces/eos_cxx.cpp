#include <eos.H>

namespace EOSData
{
  bool initialized;
  Real mintemp;
  Real maxtemp;
  Real mindens;
  Real maxdens;
  Real minx;
  Real maxx;
  Real minye;
  Real maxye;
  Real mine;
  Real maxe;
  Real minp;
  Real maxp;
  Real mins;
  Real maxs;
  Real minh;
  Real maxh;
}



// EOS initialization routine: read in general EOS parameters, then
// call any specific initialization used by the EOS.

void eos_cxx_init() {

  // Allocate and set default values

  EOSData::mintemp = 1.e-200;
  EOSData::maxtemp = 1.e200;
  EOSData::mindens = 1.e-200;
  EOSData::maxdens = 1.e200;
  EOSData::minx = 1.e-200;
  EOSData::maxx = 1.0 + 1.e-12;
  EOSData::minye = 1.e-200;
  EOSData::maxye = 1.0 + 1.e-12;
  EOSData::mine = 1.e-200;
  EOSData::maxe = 1.e200;
  EOSData::minp = 1.e-200;
  EOSData::maxp = 1.e200;
  EOSData::mins = 1.e-200;
  EOSData::maxs = 1.e200;
  EOSData::minh = 1.e-200;
  EOSData::maxh = 1.e200;

  // Set up any specific parameters or initialization steps required by the EOS we are using.
  actual_eos_cxx_init();

  // we take the approach that the Fortran initialization is the
  // reference and it has already been done, so we get any overrides
  // to these from Fortran, to ensure we are consistent.


  // If they exist, save the minimum permitted user temperature and density.
  // These are only relevant to this module if they are larger than the minimum
  // possible EOS quantities. We will reset them to be equal to the EOS minimum
  // if they are smaller than that.

  Real scratch;
  eos_get_small_temp(&scratch);
  EOSData::mintemp = scratch;

  eos_get_small_dens(&scratch);
  EOSData::mindens = scratch;

  EOSData::initialized = true

}


void eos_cxx_finalize() {

  actual_eos_cxx_finalize();

}


void eos_cxx(const eos_input_t input, eos_t& state, bool use_raw_inputs) {
  
  // Input arguments

  bool has_been_reset = false;
  bool use_composition_routine = true;


  // Local variables

#ifndef AMREX_USE_GPU
  if (!EOSData::initialized) {
    amrex::Error("EOS: not initialized");
  }
#endif

  if (use_raw_inputs) {
    use_composition_routine = false;
  }

  if (use_composition_routine) {
    // Get abar, zbar, etc.
    composition(state);
  }

  // Force the inputs to be valid.
  reset_inputs(input, state, has_been_reset);

  // Allow the user to override any details of the
  // EOS state. This should generally occur right
  // before the actual_eos call.

  //call eos_override(state)

  // Call the EOS.

  if (!has_been_reset) {
    call actual_eos_cxx(input, state);
  }
}

  subroutine reset_inputs(input, state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, &
                               eos_input_rt, eos_input_re, eos_input_rh, eos_input_tp, &
                               eos_input_rp, eos_input_th, eos_input_ph, eos_input_ps

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    // Reset the input quantities to valid values. For inputs other than rho and T,
    // this will evolve an EOS call, which will negate the need to do the main EOS call.

    if (input .eq. eos_input_rt) then

       call reset_rho(state, has_been_reset)
       call reset_T(state, has_been_reset)

    elseif (input .eq. eos_input_rh) then

       call reset_rho(state, has_been_reset)
       call reset_h(state, has_been_reset)

    elseif (input .eq. eos_input_tp) then

       call reset_T(state, has_been_reset)
       call reset_p(state, has_been_reset)

    elseif (input .eq. eos_input_rp) then

       call reset_rho(state, has_been_reset)
       call reset_p(state, has_been_reset)

    elseif (input .eq. eos_input_re) then

       call reset_rho(state, has_been_reset)
       call reset_e(state, has_been_reset)

    elseif (input .eq. eos_input_ps) then

       call reset_p(state, has_been_reset)
       call reset_s(state, has_been_reset)

    elseif (input .eq. eos_input_ph) then

       call reset_p(state, has_been_reset)
       call reset_h(state, has_been_reset)

    elseif (input .eq. eos_input_th) then

       call reset_t(state, has_been_reset)
       call reset_h(state, has_been_reset)

    endif

  end subroutine reset_inputs



  // For density, just ensure that it is within mindens and maxdens.

  subroutine reset_rho(state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, mindens, maxdens

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    state % rho = min(maxdens, max(mindens, state % rho))

  end subroutine reset_rho



  // For temperature, just ensure that it is within mintemp and maxtemp.

  subroutine reset_T(state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, mintemp, maxtemp

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    state % T = min(maxtemp, max(mintemp, state % T))

  end subroutine reset_T



  subroutine reset_e(state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, mine, maxe

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    if (state % e .lt. mine .or. state % e .gt. maxe) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_e



  subroutine reset_h(state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, minh, maxh

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    if (state % h .lt. minh .or. state % h .gt. maxh) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_h



  subroutine reset_s(state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, mins, maxs

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    if (state % s .lt. mins .or. state % s .gt. maxs) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_s



  subroutine reset_p(state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, minp, maxp

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    if (state % p .lt. minp .or. state % p .gt. maxp) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_p



  // Given an EOS state, ensure that rho and T are
  // valid, then call with eos_input_rt.

  subroutine eos_reset(state, has_been_reset)

    !$acc routine seq

    use eos_type_module, only: eos_t, eos_input_rt, mintemp, maxtemp, mindens, maxdens
    use actual_eos_module, only: actual_eos

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    !$gpu

    state % T = min(maxtemp, max(mintemp, state % T))
    state % rho = min(maxdens, max(mindens, state % rho))

    call actual_eos(eos_input_rt, state)

    has_been_reset = .true.

  end subroutine eos_reset



#ifndef AMREX_USE_GPU
  subroutine check_inputs(input, state)

    use amrex_error_module
    use network, only: nspec
    use eos_type_module, only: eos_t, print_state, minx, maxx, minye, maxye, &
                               eos_input_rt, eos_input_re, eos_input_rp, eos_input_rh, &
                               eos_input_th, eos_input_tp, eos_input_ph, eos_input_ps

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    integer :: n

    // Check the inputs for validity.

    do n = 1, nspec
       if (state % xn(n) .lt. minx) then
          call print_state(state)
          call amrex_error('EOS: mass fraction less than minimum possible mass fraction.')
       else if (state % xn(n) .gt. maxx) then
          call print_state(state)
          call amrex_error('EOS: mass fraction more than maximum possible mass fraction.')
       endif
    enddo

    if (state % y_e .lt. minye) then
       call print_state(state)
       call amrex_error('EOS: y_e less than minimum possible electron fraction.')
    else if (state % y_e .gt. maxye) then
       call print_state(state)
       call amrex_error('EOS: y_e greater than maximum possible electron fraction.')
    endif

    if (input .eq. eos_input_rt) then

       call check_rho(state)
       call check_T(state)

    elseif (input .eq. eos_input_rh) then

       call check_rho(state)
       call check_h(state)

    elseif (input .eq. eos_input_tp) then

       call check_T(state)
       call check_p(state)

    elseif (input .eq. eos_input_rp) then

       call check_rho(state)
       call check_p(state)

    elseif (input .eq. eos_input_re) then

       call check_rho(state)
       call check_e(state)

    elseif (input .eq. eos_input_ps) then

       call check_p(state)
       call check_s(state)

    elseif (input .eq. eos_input_ph) then

       call check_p(state)
       call check_h(state)

    elseif (input .eq. eos_input_th) then

       call check_t(state)
       call check_h(state)

    endif

  end subroutine check_inputs



  subroutine check_rho(state)

    use amrex_error_module
    use eos_type_module, only: eos_t, mindens, maxdens, print_state

    implicit none

    type (eos_t), intent(in) :: state

    if (state % rho .lt. mindens) then
       call print_state(state)
       call amrex_error('EOS: rho smaller than mindens.')
    else if (state % rho .gt. maxdens) then
       call print_state(state)
       call amrex_error('EOS: rho greater than maxdens.')
    endif

  end subroutine check_rho



  subroutine check_T(state)

    use amrex_error_module
    use eos_type_module, only: eos_t, mintemp, maxtemp, print_state

    implicit none

    type (eos_t), intent(in) :: state

    if (state % T .lt. mintemp) then
       call print_state(state)
       call amrex_error('EOS: T smaller than mintemp.')
    else if (state % T .gt. maxtemp) then
       call print_state(state)
       call amrex_error('EOS: T greater than maxtemp.')
    endif

  end subroutine check_T



  subroutine check_e(state)

    use amrex_error_module
    use eos_type_module, only: eos_t, mine, maxe, print_state

    implicit none

    type (eos_t), intent(in) :: state

    if (state % e .lt. mine) then
       call print_state(state)
       call amrex_error('EOS: e smaller than mine.')
    else if (state % e .gt. maxe) then
       call print_state(state)
       call amrex_error('EOS: e greater than maxe.')
    endif

  end subroutine check_e



  subroutine check_h(state)

    use amrex_error_module
    use eos_type_module, only: eos_t, minh, maxh, print_state

    implicit none

    type (eos_t), intent(in) :: state

    if (state % h .lt. minh) then
       call print_state(state)
       call amrex_error('EOS: h smaller than minh.')
    else if (state % h .gt. maxh) then
       call print_state(state)
       call amrex_error('EOS: h greater than maxh.')
    endif

  end subroutine check_h



  subroutine check_s(state)

    use amrex_error_module
    use eos_type_module, only: eos_t, mins, maxs, print_state

    implicit none

    type (eos_t), intent(in) :: state

    if (state % s .lt. mins) then
       call print_state(state)
       call amrex_error('EOS: s smaller than mins.')
    else if (state % s .gt. maxs) then
       call print_state(state)
       call amrex_error('EOS: s greater than maxs.')
    endif

  end subroutine check_s



  subroutine check_p(state)

    use amrex_error_module
    use eos_type_module, only: eos_t, minp, maxp, print_state

    implicit none

    type (eos_t), intent(in) :: state

    if (state % p .lt. minp) then
       call print_state(state)
       call amrex_error('EOS: p smaller than minp.')
    else if (state % p .gt. maxp) then
       call print_state(state)
       call amrex_error('EOS: p greater than maxp.')
    endif

  end subroutine check_p
#endif

end module eos_module
