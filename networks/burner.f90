module burner_module

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use specific_burner_module

  interface burner
     module procedure scalar_burner
     module procedure vector_burner
     module procedure rank2_burner
     module procedure rank3_burner
  end interface burner
  
contains

  subroutine scalar_burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t), intent(inout) :: state_in
    type (eos_t), intent(inout) :: state_out
    double precision, intent(in) :: dt, time

    type (eos_t) :: vector_state_in(1), vector_state_out(1)

    call vector_burner(vector_state_in, vector_state_out, dt, time)

    state_in  = vector_state_in(1)
    state_out = vector_state_out(1)
    
  end subroutine scalar_burner


  
  subroutine vector_burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t), intent(inout)  :: state_in(:)
    type (eos_t), intent(inout)  :: state_out(:)
    double precision, intent(in) :: dt, time

    integer :: j, N

    N = size(state_in)
    
    ! Make sure the network has been initialized.
    
    if (.NOT. network_initialized) then
       call bl_error("ERROR in burner: must initialize network first.")
    endif

    ! We assume that the valid quantities coming in are (rho, e); do an EOS call
    ! to make sure all other variables are consistent.

    call eos(eos_input_re, state_in)

    ! Initialize the final state by assuming it does not change.

    state_out = state_in

    call specific_burner(state_in, state_out, dt, time)

    do j = 1, N

       ! Store the new mass fractions -- note, we discard the temperature
       ! here.  Make sure that they are positive and less than one.
       
       state_out(j) % xn(:) = max(smallx, min(ONE, state_out(j) % xn(:)))

       ! Enforce sum{X_k} = 1.

       state_out(j) % xn(:) = state_out(j) % xn(:) / sum(state_out(j) % xn(:))

    enddo
       
    ! Now update the temperature to match the new internal energy.

    call eos(eos_input_re, state_out)

  end subroutine vector_burner


  
  subroutine rank2_burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t), intent(inout) :: state_in(:,:)
    type (eos_t), intent(inout) :: state_out(:,:)
    double precision, intent(in) :: dt, time

    type (eos_t) :: vector_state_in(size(state_in)), vector_state_out(size(state_out))

    vector_state_in(:)  = reshape(state_in(:,:), (/ size(state_in) /))
    vector_state_out(:) = reshape(state_out(:,:), (/ size(state_out) /))
    
    call vector_burner(vector_state_in, vector_state_out, dt, time)

    state_in(:,:)  = reshape(vector_state_in, (/ size(state_in,1),size(state_in,2) /))
    state_out(:,:) = reshape(vector_state_out, (/ size(state_out,1),size(state_out,2) /))
    
  end subroutine rank2_burner



  subroutine rank3_burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t), intent(inout) :: state_in(:,:,:)
    type (eos_t), intent(inout) :: state_out(:,:,:)
    double precision, intent(in) :: dt, time

    type (eos_t) :: vector_state_in(size(state_in)), vector_state_out(size(state_out))

    vector_state_in(:)  = reshape(state_in(:,:,:), (/ size(state_in) /))
    vector_state_out(:) = reshape(state_out(:,:,:), (/ size(state_out) /))
    
    call vector_burner(vector_state_in, vector_state_out, dt, time)

    state_in(:,:,:)  = reshape(vector_state_in, (/ size(state_in,1),size(state_in,2),size(state_in,3) /))
    state_out(:,:,:) = reshape(vector_state_out, (/ size(state_out,1),size(state_out,2),size(state_out,3) /))
    
  end subroutine rank3_burner

end module burner_module
