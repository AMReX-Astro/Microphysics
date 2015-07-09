module eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module
  use specific_eos_module

  implicit none

  public eos_init, eos

  interface eos
     module procedure scalar_eos
     module procedure vector_eos
     module procedure rank2_eos
     module procedure rank3_eos
     module procedure rank4_eos
     module procedure rank5_eos
  end interface eos

contains

  ! EOS initialization routine: read in general EOS parameters, then 
  ! call any specific initialization used by the EOS.

  subroutine eos_init(small_temp, small_dens)

    implicit none
 
    double precision, intent(in), optional :: small_temp
    double precision, intent(in), optional :: small_dens
 
    smallt = 1.d4

    if (present(small_temp)) then
      if (small_temp > ZERO) then
       smallt = small_temp
      end if
    endif

    smalld = 1.d-5
 
    if (present(small_dens)) then
       if (small_dens > ZERO) then
         smalld = small_dens
       endif
    endif

    ! Set up any specific parameters or initialization steps required by the EOS we are using.

    call specific_eos_init

  end subroutine eos_init



  subroutine scalar_eos(input, state, do_eos_diag)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state
    logical, optional, intent(in   ) :: do_eos_diag

    type (eos_t) :: vector_state(1)

    vector_state(1) = state

    if (present(do_eos_diag)) then

       call vector_eos(input, vector_state, do_eos_diag)

    else

       call vector_eos(input, vector_state)

    endif

    state = vector_state(1)

  end subroutine scalar_eos


  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------

  subroutine vector_eos(input, state, do_eos_diag)

    implicit none

    ! Input arguments

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state(:)
    logical, optional, intent(in   ) :: do_eos_diag

    integer :: j, N

    ! Local variables and arrays

    double precision :: ymass(nspec), ysum, yzsum
    double precision :: e_want, p_want, s_want, h_want

    logical eosfail, eos_diag

    integer :: ns, ierr

    N = size(state)
    
    if (.not. initialized) call bl_error('EOS: not initialized')

    eos_diag = .false.

    if (present(do_eos_diag)) eos_diag = do_eos_diag

    do j = 1, N
       ! Check to make sure the composition was set properly.

       do ns = 1, nspec
         if (state(j) % xn(ns) .lt. init_test) call eos_error(ierr_init_xn, input, state(j) % loc)
       enddo

       ! Get abar, zbar, etc.

       call composition(state(j), .false.)
    enddo

    eosfail = .false.

    ierr = 0

    ! Check the inputs, and do initial setup for iterations

    do j = 1, N

       if ( state(j) % T .lt. mintemp ) then
          print *, 'TEMP = ', state(j) % T
          call bl_error('EOS: temp less than minimum possible temperature.')
          eosfail = .true.
          return
       end if
       if ( state(j) % T .gt. maxtemp ) then
          print *, 'TEMP = ', state(j) % T
          call bl_error('EOS: temp greater than maximum possible temperature.')
          eosfail = .true.
          return
       end if

       if ( state(j) % rho .lt. mindens ) then
          print *, 'DENS = ', state(j) % rho
          call bl_error('EOS: dens less than minimum possible density.')
          eosfail = .true.
          return
       end if
       if ( state(j) % rho .gt. maxdens ) then
          print *, 'DENS = ', state(j) % rho
          call bl_error('EOS: dens greater than maximum possible density.')
          eosfail = .true.
       end if

       if (input .eq. eos_input_rt) then

         if (state(j) % rho .lt. init_test .or. state(j) % T .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       elseif (input .eq. eos_input_rh) then

         if (state(j) % rho .lt. init_test .or. state(j) % h .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       elseif (input .eq. eos_input_tp) then

         if (state(j) % T   .lt. init_test .or. state(j) % p .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       elseif (input .eq. eos_input_rp) then

         if (state(j) % rho .lt. init_test .or. state(j) % p .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       elseif (input .eq. eos_input_re) then

         if (state(j) % rho .lt. init_test .or. state(j) % e .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       elseif (input .eq. eos_input_ps) then

         if (state(j) % p   .lt. init_test .or. state(j) % s .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       elseif (input .eq. eos_input_ph) then

         if (state(j) % p   .lt. init_test .or. state(j) % h .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       elseif (input .eq. eos_input_th) then

         if (state(j) % T   .lt. init_test .or. state(j) % h .lt. init_test) call eos_error(ierr_init, input, state(j) % loc)

       endif

    enddo

    ! Call the EOS.

    call specific_eos(eosfail, state, N, input)

    ! Get dpdX, dedX, dhdX.

    do j = 1, N
       call composition_derivatives(state(j), .false.)
    enddo

  end subroutine vector_eos



  subroutine rank2_eos(input, state, do_eos_diag)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state(:,:)
    logical, optional, intent(in   ) :: do_eos_diag

    type (eos_t) :: vector_state(size(state))

    vector_state(:) = reshape(state(:,:), (/ size(state) /))

    if (present(do_eos_diag)) then

       call vector_eos(input, vector_state, do_eos_diag)

    else

       call vector_eos(input, vector_state)

    endif

    state(:,:) = reshape(vector_state, (/ size(state,1),size(state,2) /))

  end subroutine rank2_eos



  subroutine rank3_eos(input, state, do_eos_diag)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state(:,:,:)
    logical, optional, intent(in   ) :: do_eos_diag

    type (eos_t) :: vector_state(size(state))

    vector_state(:) = reshape(state(:,:,:), (/ size(state) /))

    if (present(do_eos_diag)) then

       call vector_eos(input, vector_state, do_eos_diag)

    else

       call vector_eos(input, vector_state)

    endif

    state(:,:,:) = reshape(vector_state, (/ size(state,1),size(state,2),size(state,3) /))

  end subroutine rank3_eos



  subroutine rank4_eos(input, state, do_eos_diag)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state(:,:,:,:)
    logical, optional, intent(in   ) :: do_eos_diag

    type (eos_t) :: vector_state(size(state))

    vector_state(:) = reshape(state(:,:,:,:), (/ size(state) /))

    if (present(do_eos_diag)) then

       call vector_eos(input, vector_state, do_eos_diag)

    else

       call vector_eos(input, vector_state)

    endif

    state(:,:,:,:) = reshape(vector_state, (/ size(state,1),size(state,2),size(state,3),size(state,4) /))

  end subroutine rank4_eos



  subroutine rank5_eos(input, state, do_eos_diag)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state(:,:,:,:,:)
    logical, optional, intent(in   ) :: do_eos_diag

    type (eos_t) :: vector_state(size(state))

    vector_state(:) = reshape(state(:,:,:,:,:), (/ size(state) /))

    if (present(do_eos_diag)) then

       call vector_eos(input, vector_state, do_eos_diag)

    else

       call vector_eos(input, vector_state)

    endif

    state(:,:,:,:,:) = reshape(vector_state, (/ size(state,1),size(state,2),size(state,3),size(state,4),size(state,5) /))

  end subroutine rank5_eos



  subroutine eos_error(err, input, pt_index)

    implicit none

    integer,           intent(in) :: err
    integer,           intent(in) :: input
    integer, optional, intent(in) :: pt_index(3)

    integer :: dim_ptindex

    character (len=64) :: err_string, zone_string, eos_input_str

    write(eos_input_str, '(A13, I1)') ' EOS input = ', input

    if (err .eq. ierr_general) then

      err_string = 'EOS: error in the EOS.'

    elseif (err .eq. ierr_input) then

      err_string = 'EOS: invalid input.'

    elseif (err .eq. ierr_iter_conv) then

      err_string = 'EOS: Newton-Raphson iterations failed to converge.'

    elseif (err .eq. ierr_neg_e) then

      err_string = 'EOS: energy < 0 in the EOS.'

    elseif (err .eq. ierr_neg_p) then

      err_string = 'EOS: pressure < 0 in the EOS.'

    elseif (err .eq. ierr_neg_h) then

      err_string = 'EOS: enthalpy < 0 in the EOS.'

    elseif (err .eq. ierr_neg_s) then

      err_string = 'EOS: entropy < 0 in the EOS.'

    elseif (err .eq. ierr_init) then

      err_string = 'EOS: the input variables were not initialized.'

    elseif (err .eq. ierr_init_xn) then

      err_string = 'EOS: the species abundances were not initialized.'

    elseif (err .eq. ierr_iter_var) then

      err_string = 'EOS: the variable you are iterating over was not recognized.'

    else

      err_string = 'EOS: invalid input to error handler.'

    endif

    err_string = err_string // eos_input_str

    ! this format statement is for writing into zone_string -- make sure that
    ! the len of z_err can accomodate this format specifier
1001 format(1x,"zone index info: i = ", i5)
1002 format(1x,"zone index info: i = ", i5, '  j = ', i5)
1003 format(1x,"zone index info: i = ", i5, '  j = ', i5, '  k = ', i5)

    if (present(pt_index)) then

       dim_ptindex = 3

       if (pt_index(3) .eq. -99) then
          dim_ptindex = 2
          if (pt_index(2) .eq. -99) then
             dim_ptindex = 1
             if (pt_index(1) .eq. -99) then
                dim_ptindex = 0
             endif
          endif
       endif

       if (dim_ptindex .eq. 1) then 
          write (zone_string,1001) pt_index(1)
       else if (dim_ptindex .eq. 2) then 
          write (zone_string,1002) pt_index(1), pt_index(2)
       else if (dim_ptindex .eq. 3) then 
          write (zone_string,1003) pt_index(1), pt_index(2), pt_index(3)
       end if

    else

      zone_string = ''

    endif

    call bl_error(err_string, zone_string)

  end subroutine eos_error


end module eos_module
