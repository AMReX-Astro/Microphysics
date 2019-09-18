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


  subroutine nonaka_rhs(state)

    use extern_probin_module, only: nonaka_i, nonaka_j, nonaka_k, nonaka_file
    use amrex_fort_module, only: rt => amrex_real
    use burn_type_module, only: burn_t
    use actual_network, only: nspec, aion

    implicit none

    type (burn_t), intent(in) :: state

    integer :: nonaka_file_unit, j

    character(len=20) :: vector_format = ''
    character(len=20) :: scalar_format = ''
    
    if (state % i == nonaka_i .and. &
        state % j == nonaka_j .and. &
        state % k == nonaka_k) then

        ! append current state to nonaka log
        ! at the current simulation time

        write(vector_format, '("(", I0, "E30.16E5", ")")') nspec
        write(scalar_format, '("(", I0, "E30.16E5", ")")') 1
        
        open(newunit=nonaka_file_unit, file=nonaka_file, status="old", position="append", action="write")
        write(unit=nonaka_file_unit, fmt=scalar_format, advance="no") state % time

        ! Mass fractions X
        write(unit=nonaka_file_unit, fmt=vector_format, advance="no") (state % xn(j), j = 1, nspec)

        ! Convert molar fraction rhs to mass fraction rhs dX/dt
        write(unit=nonaka_file_unit, fmt=vector_format) (state % ydot(j) * aion(j), j = 1, nspec)
        close(unit=nonaka_file_unit)

    end if
      
  end subroutine nonaka_rhs

end module nonaka_plot_module
