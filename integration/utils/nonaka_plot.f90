module nonaka_plot_module

  implicit none

  public

contains

  subroutine nonaka_init()

    use extern_probin_module, only: nonaka_file
    use actual_network, only: nspec, short_spec_names

    implicit none

    integer :: nonaka_file_unit, i
    character(len=20*nspec+10) :: header_line
    character(len=20) :: scratch

    header_line = "time"

    do i = 1, nspec
      scratch = "X("//trim(short_spec_names(i))//")"
      header_line = trim(header_line)//" "//trim(scratch)
    end do

    do i = 1, nspec
      scratch = "dXdt("//trim(short_spec_names(i))//")"
      header_line = trim(header_line)//" "//trim(scratch)
    end do

    open(newunit=nonaka_file_unit, file=nonaka_file, status="replace", action="write")
    write(unit=nonaka_file_unit, fmt=*) trim(header_line)
    close(unit=nonaka_file_unit)

  end subroutine nonaka_init


  subroutine nonaka_rhs(state, reference_time, trim_after_timestep)

    ! state: the burn_t corresponding to the current state
    !        with state % time relative to the start of the current burn call.
    ! reference_time: the simulation time at the start of the current burn call.
    !
    ! The current simulation time is state % time + reference_time
    !
    ! trim_after_timestep: if trim_after_timestep = .true. then this call is after the
    !               ODE integrator has finished the timestep for the current burn call.
    !               In that case, we trim entries in the nonaka file past the timestep end.

    use extern_probin_module, only: nonaka_i, nonaka_j, nonaka_k, nonaka_file
    use amrex_fort_module, only: rt => amrex_real
    use burn_type_module, only: burn_t
    use actual_network, only: nspec, aion

    implicit none

    type (burn_t), intent(in) :: state
    real(rt),      intent(in) :: reference_time
    
    ! optional: trim entries past end of timestep
    logical, intent(in), optional :: trim_after_timestep

    integer :: nonaka_file_unit, j
    integer :: i, nextline
    real(rt) :: tprev

    character(len=20) :: vector_format = ''
    character(len=20) :: scalar_format = ''

    if (state % i == nonaka_i .and. &
        state % j == nonaka_j .and. &
        state % k == nonaka_k) then

        ! append current state to nonaka log
        ! at the current simulation time

        write(vector_format, '("(", I0, "E30.16E5", ")")') nspec
        write(scalar_format, '("(", I0, "E30.16E5", ")")') 1
        
        open(newunit=nonaka_file_unit, file=nonaka_file, status="old", position="append", action="readwrite", &
             access="stream", form="formatted")

        inquire(unit=nonaka_file_unit, pos=nextline)
        
        if (present(trim_after_timestep)) then
            if (trim_after_timestep) then
                ! determine last timestep in file
                nextline = nextline - ( 30*(1 + 2*nspec) + 1 )
                read(unit=nonaka_file_unit, pos=nextline, fmt=scalar_format) tprev

                i = 0
                ! remove last rows where VODE took a much too large timestep
                do while (tprev >= time .and. i < 2)
                    nextline = nextline - ( 30*(1 + 2*nspec) + 1 )
                    read(unit=nonaka_file_unit, pos=nextline, fmt=scalar_format) tprev
                    i = i + 1
                end do

                ! store where to write new data
                nextline = nextline + ( 30*(1 + 2*nspec) + 1 )
            end if
        end if
           
        write(unit=nonaka_file_unit, fmt=scalar_format, pos=nextline, advance="no") (state % time + reference_time)

        ! Mass fractions X
        write(unit=nonaka_file_unit, fmt=vector_format, advance="no") (state % xn(j), j = 1, nspec)

        ! Convert molar fraction rhs to mass fraction rhs dX/dt
        write(unit=nonaka_file_unit, fmt=vector_format) (state % ydot(j) * aion(j), j = 1, nspec)
        close(unit=nonaka_file_unit)

    end if

  end subroutine nonaka_rhs

end module nonaka_plot_module
