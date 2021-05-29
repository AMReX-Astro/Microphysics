module actual_network

  use network_properties

  integer, parameter :: nrates = 1
  integer, parameter :: num_rate_groups = 1

  character (len=32), parameter :: network_name = "vode_example"

contains

  subroutine actual_network_init

    implicit none

    call network_properties_init()

  end subroutine actual_network_init


  subroutine actual_network_finalize

    implicit none

    call network_properties_finalize()

  end subroutine actual_network_finalize

end module actual_network
