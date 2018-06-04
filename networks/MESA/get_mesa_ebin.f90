!
! routine to return the binding energy using MESA
!

function get_mesa_ebin(i) result(ebin)

   ! BoxLib
   use bl_types,  only: dp_t
   use bl_error_module, only: bl_error

   ! MAESTRO
   use network,   only: chem_id

   ! MESA
   use net_lib,   only: chem_isos, del_Mn, del_Mp
   use const_def, only: Qconv

   implicit none

   ! INPUT
   integer, intent(in) :: i

   ! OUTPUT
   real(kind=dp_t) :: ebin

   ! LOCAL
   integer :: cid

   cid = chem_id(i)

   ! Qconv converts from MeV to erg and multiplies by N_A = 6.024d23
   ! Dividing by atomic weight results in [ebin] = erg/g
   !ebin = chem_isos%binding_energy(cid) / chem_isos%W(cid)
   ebin = (chem_isos%binding_energy(cid) - chem_isos%Z(cid)*del_Mp - &
            chem_isos%N(cid)*del_Mn) / chem_isos%Z_plus_N(cid)

   ebin = ebin*Qconv

   return

end function get_mesa_ebin

