Module xnet_interface
  Implicit None
  
  Interface benuc
    Subroutine benuc(y,enb,enm,ytot,ztot,atot)
      Use xnet_types, Only: dp
      Real(dp), Intent(in) :: y(:)
      Real(dp), Intent(out) :: enb
      Real(dp), Intent(out) :: enm
      Real(dp), Intent(out) :: ytot, ztot, atot
    End Subroutine benuc
  End Interface

  Interface cross_sect
    Subroutine cross_sect(mask_in)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine cross_sect
  End Interface

  Interface ffn_rate
    Subroutine ffn_rate(nffn,t9,ene,rf,dlnrfdt9,mask_in)
      Use xnet_types, Only: dp
      Integer, Intent(in) :: nffn
      Real(dp), Intent(in) :: t9(:)
      Real(dp), Intent(in) :: ene(:)
      Real(dp), Intent(out) :: rf(nffn,size(t9))
      Real(dp), Intent(out) :: dlnrfdt9(nffn,size(t9))
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine ffn_rate
  End Interface

  Interface final_output
    Subroutine final_output(kstep,mask_in)
      Integer, Intent(in) :: kstep
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine final_output
  End Interface

  Interface flux
    Subroutine flux(mask_in)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine flux
  End Interface

  Interface flux_check
    Subroutine flux_check(mask_in)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine flux_check
  End Interface

  Interface flux_init
    Subroutine flux_init
    End Subroutine flux_init
  End Interface
  
  Interface flux_search
    Subroutine flux_search(iflx,iw,i7,q,desc,mflx)
      Use xnet_types, Only: dp
      Integer, Intent(in) :: iw
      Integer, Intent(in) :: i7(7)
      Real(dp), Intent(in) :: q
      Character(4), Intent(in) :: desc
      Integer, Intent(inout) :: mflx
      Integer, Intent(out) :: iflx
    End Subroutine flux_search
  End Interface

  Interface format_nuclear_data
    Subroutine format_nuclear_data(lun_out,data_dir)
      Integer, Intent(in) :: lun_out
      Character(*), Intent(in) :: data_dir
    End Subroutine format_nuclear_data
  End Interface

  Interface full_net
    Subroutine full_net(kstep)
      Integer, Intent(out) :: kstep
    End Subroutine full_net
  End Interface

  Interface index_from_name
    Subroutine index_from_name(nuc_name,inuc)
      Character(*), Intent(in) :: nuc_name
      Integer, Intent(out) :: inuc
    End Subroutine index_from_name
  End Interface

  Interface jacobian_bksub
    Subroutine jacobian_bksub(kstep,yrhs,dy,t9rhs,dt9,mask_in)
      Use xnet_types, Only: dp
      Integer, Intent(in) :: kstep
      Real(dp), Intent(in) :: yrhs(:,:)
      Real(dp), Intent(in) :: t9rhs(:)
      Real(dp), Intent(out) :: dy(size(yrhs,1),size(yrhs,2))
      Real(dp), Intent(out) :: dt9(size(t9rhs))
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine jacobian_bksub
  End Interface

  Interface jacobian_build
    Subroutine jacobian_build(diag,mult,mask_in)
      Use xnet_types, Only: dp
      Real(dp), Intent(in) :: diag(:), mult
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine jacobian_build
  End Interface

  Interface jacobian_decomp
    Subroutine jacobian_decomp(kstep,mask_in)
      Integer, Intent(in) :: kstep
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine jacobian_decomp
  End Interface

  Interface jacobian_finalize
    Subroutine jacobian_finalize
    End Subroutine jacobian_finalize
  End Interface

  Interface jacobian_solve
    Subroutine jacobian_solve(kstep,yrhs,dy,t9rhs,dt9,mask_in)
      Use xnet_types, Only: dp
      Integer, Intent(in) :: kstep
      Real(dp), Intent(in) :: yrhs(:,:)
      Real(dp), Intent(in) :: t9rhs(:)
      Real(dp), Intent(out) :: dy(size(yrhs,1),size(yrhs,2))
      Real(dp), Intent(out) :: dt9(size(t9rhs))
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine jacobian_solve
  End Interface

  Interface match_react
    Subroutine match_react(data_dir)
      Character(*), Intent(in) :: data_dir
    End Subroutine match_react
  End Interface

  Interface name_ordered
    Subroutine name_ordered(base_string,n,nmax)
      Integer, Intent(in) :: n, nmax
      Character(*), Intent(inout) :: base_string
    End Subroutine name_ordered
  End Interface

  Interface net_preprocess
    Subroutine net_preprocess(lun_out,data_dir,data_desc)
      Integer, Intent(in) :: lun_out
      Character(*), Intent(in) :: data_dir
      Character(*), Intent(in) :: data_desc
    End Subroutine net_preprocess
  End Interface

  Interface norm
    Subroutine norm(yy)
      Use xnet_types, Only: dp
      Real(dp), Intent(inout) :: yy(:)
    End Subroutine norm
  End Interface

  Interface partf
    Subroutine partf(t9,mask_in)
      Use xnet_types, Only: dp
      Real(dp), Intent(in) :: t9(:)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine partf
  End Interface

  Interface read_jacobian_data
    Subroutine read_jacobian_data(data_dir)
      Character(*), Intent(in) :: data_dir
    End Subroutine read_jacobian_data
  End Interface

  Interface read_match_data
    Subroutine read_match_data(data_dir)
      Character(*), Intent(in) :: data_dir
    End Subroutine read_match_data
  End Interface

  Interface read_nuclear_data
    Subroutine read_nuclear_data(data_dir,data_desc)
      Character(*), Intent(in) :: data_dir
      Character(80), Intent(out) :: data_desc
    End Subroutine read_nuclear_data
  End Interface

  Interface read_reaction_data
    Subroutine read_reaction_data(data_dir)
      Character(*), Intent(in) :: data_dir
    End Subroutine read_reaction_data
  End Interface

  Interface read_ffn_data
    Subroutine read_ffn_data(nffn,data_dir)
      Integer, Intent(in)  :: nffn
      Character(*), Intent(in) :: data_dir
    End Subroutine read_ffn_data
  End Interface

  Interface screening
    Subroutine screening(mask_in)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine screening
  End Interface

  Interface solve_be
    Subroutine solve_be(kstep,its)
      Integer, Intent(in) :: kstep
      Integer, Intent(inout) :: its(:)
    End Subroutine solve_be
  End Interface

  Interface sparse_check
    Subroutine sparse_check(data_dir)
      Character(*), Intent(in) :: data_dir
    End Subroutine sparse_check
  End Interface

  Interface step_be
    Subroutine step_be(kstep,inr)
      Integer, Intent(in) :: kstep
      Integer, Intent(inout) :: inr(:)
    End Subroutine step_be
  End Interface

  Interface string_lc
    Subroutine string_lc(string)
      Character(*), Intent(inout) :: string
    End Subroutine string_lc
  End Interface

  Interface t9rhofind
    Subroutine t9rhofind(kstep,tf,nf,t9f,rhof,mask_in)
      Use xnet_types, Only: dp
      Integer, Intent(in) :: kstep
      Real(dp), Intent(in) :: tf(:)
      Integer, Intent(inout) :: nf(size(tf))
      Real(dp), Intent(inout) :: t9f(size(tf)), rhof(size(tf))
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine t9rhofind
  End Interface

  Interface timestep
    Subroutine timestep(kstep,mask_in)
      Integer, Intent(in) :: kstep
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine timestep
  End Interface

  Interface ts_output
    Subroutine ts_output(kstep,enuc,edot,mask_in)
      Use xnet_types, Only: dp
      Integer, Intent(in) :: kstep
      Real(dp), Intent(in) :: enuc(:), edot(:)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine ts_output
  End Interface

  Interface update_iweak
    Subroutine update_iweak(t9,mask_in)
      Use xnet_types, Only: dp
      Real(dp), Intent(in) :: t9(:)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine update_iweak
  End Interface

  Interface xnet_terminate
    Subroutine xnet_terminate(c_diagnostic,i_diagnostic)
      Character(*), Intent(in) :: c_diagnostic
      Integer, Intent(in), Optional :: i_diagnostic
    End Subroutine xnet_terminate
  End Interface

  Interface yderiv
    Subroutine yderiv(mask_in)
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine yderiv
  End Interface

  Interface yderiv_and_jacobian_build
    Subroutine yderiv_and_jacobian_build(diag,mult,mask_in)
      Use xnet_types, Only: dp
      Real(dp), Intent(in) :: diag(:), mult
      Logical, Optional, Target, Intent(in) :: mask_in(:)
    End Subroutine yderiv_and_jacobian_build
  End Interface

  Interface ye_norm
    Subroutine ye_norm(yy,ye)
      Use xnet_types, Only: dp
      Real(dp), Intent(in) :: ye
      Real(dp), Intent(inout) :: yy(:)
    End Subroutine ye_norm
  End Interface

  Interface
    Subroutine xnet_init(data_dir,data_desc)
      Character(*), Intent(in) :: data_dir
      Character(80), Intent(out) :: data_desc
    End Subroutine bn_xnetInit
  End Interface

  Interface
    Subroutine xnet_finalize()
    End Subroutine xnet_finalize
  End Interface

End Module xnet_interface