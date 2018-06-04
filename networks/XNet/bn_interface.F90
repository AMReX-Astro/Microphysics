Module bn_interface

  interface
     subroutine bn_xnetInit(data_dir,data_desc)
        implicit none
        character(*), intent(in) :: data_dir
        character(80), intent(out) :: data_desc
     end subroutine bn_xnetInit
  end interface

  interface
     subroutine bn_xnetFinalize()
        implicit none
     end subroutine bn_xnetFinalize
  end interface

  interface
     subroutine bn_burner(tstep,temp,density,xIn,xOut,sdotRate,burnedZone,kstep)
       implicit none
       logical, intent(IN), dimension(:)              :: burnedZone
       real, intent(IN)                               :: tstep
       real, intent(IN), dimension(:)                 :: temp,density
       real, intent(OUT), dimension(size(temp))       :: sdotRate
       real, intent(IN), dimension(:,:)               :: xIn
       real, intent(OUT), dimension(size(xIn,1),size(xIn,2)) :: xOut
       integer, intent(OUT)                           :: kstep
     end subroutine bn_burner
  end interface

end Module bn_interface
