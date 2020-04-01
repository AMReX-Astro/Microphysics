subroutine init_unit_test(name, namlen) bind(C, name="init_unit_test")

  use amrex_fort_module, only: rt => amrex_real
  use extern_probin_module
  use microphysics_module
  use screening_module, only: add_screening_factor

  implicit none

  integer, intent(in) :: namlen
  integer, intent(in) :: name(namlen)

  call runtime_init(name, namlen)

  call microphysics_init(small_temp, small_dens)

  ! initialize the screening factors
  call add_screening_factor(zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ihe4),aion(ihe4),4.0e0_rt,8.0e0_rt)

  call add_screening_factor(zion(ic12),aion(ic12),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ic12),aion(ic12),zion(ic12),aion(ic12))

  call add_screening_factor(zion(ic12),aion(ic12),zion(io16),aion(io16))

  call add_screening_factor(zion(io16),aion(io16),zion(io16),aion(io16))

  call add_screening_factor(zion(io16),aion(io16),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ine20),aion(ine20),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(img24),aion(img24),zion(ihe4),aion(ihe4))

  call add_screening_factor(13.0e0_rt,27.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(isi28),aion(isi28),zion(ihe4),aion(ihe4))

  call add_screening_factor(15.0e0_rt,31.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(is32),aion(is32),zion(ihe4),aion(ihe4))

  call add_screening_factor(17.0e0_rt,35.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(iar36),aion(iar36),zion(ihe4),aion(ihe4))

  call add_screening_factor(19.0e0_rt,39.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ica40),aion(ica40),zion(ihe4),aion(ihe4))

  call add_screening_factor(21.0e0_rt,43.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(iti44),aion(iti44),zion(ihe4),aion(ihe4))

  call add_screening_factor(23.0e0_rt,47.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(icr48),aion(icr48),zion(ihe4),aion(ihe4))

  call add_screening_factor(25.0e0_rt,51.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ife52),aion(ife52),zion(ihe4),aion(ihe4))

  call add_screening_factor(27.0e0_rt,55.0e0_rt,1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ife54),aion(ife54),1.0e0_rt,1.0e0_rt)

  call add_screening_factor(zion(ife54),aion(ife54),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ife56),aion(ife56),1.0e0_rt,1.0e0_rt)

  call add_screening_factor(1.0e0_rt,2.0e0_rt,zion(ih1),aion(ih1))

  call add_screening_factor(zion(ih1),aion(ih1),zion(ih1),aion(ih1))

  call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3))

  call add_screening_factor(zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4))

  call add_screening_factor(zion(ic12),aion(ic12),zion(ih1),aion(ih1))

  call add_screening_factor(zion(in14),aion(in14),zion(ih1),aion(ih1))

  call add_screening_factor(zion(io16),aion(io16),zion(ih1),aion(ih1))

  call add_screening_factor(zion(in14),aion(in14),zion(ihe4),aion(ihe4))

end subroutine init_unit_test
