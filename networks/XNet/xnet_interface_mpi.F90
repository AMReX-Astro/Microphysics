Module xnet_interface_mpi
  Implicit None

  Interface control_bcast
    Subroutine control_bcast(data_dir)
      Character(80), Intent(out) :: data_dir
    End Subroutine control_bcast
  End Interface

  Interface jacobian_bcast
    Subroutine jacobian_bcast(data_dir)
      Character(*), Intent(in) :: data_dir
    End Subroutine jacobian_bcast
  End Interface

  Interface netdata_bcast
    Subroutine netdata_bcast(data_dir,data_desc)
      Character(*), Intent(in) :: data_dir
      Character(80), Intent(out) :: data_desc
    End Subroutine netdata_bcast
  End Interface

  Interface match_bcast
    Subroutine match_bcast(data_dir)
      Character(*), Intent(in) :: data_dir
    End Subroutine match_bcast
  End Interface

End Module xnet_interface_mpi