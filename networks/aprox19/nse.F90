module nse_module

  use amrex_fort_module, only : rt => amrex_real

  use network_properties, only : nspec
  implicit none

  integer, parameter :: ntemp = 71
  integer, parameter :: nden = 31
  integer, parameter :: nye = 21

  integer, parameter :: npts = 46221


  real(rt), allocatable :: ttlog(:), ddlog(:), yetab(:)
  real(rt), allocatable :: helium(:), sica(:), fegroup(:)
  real(rt), allocatable :: abartab(:), ebtab(:), wratetab(:)
  real(rt), allocatable :: massfractab(:, :)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: ttlog, ddlog, yetab
  attributes(managed) :: helium, sica, fegroup
  attributes(managed) :: abartab, ebtab, wratetab
  attributes(managed) :: massfractab
#endif

  !$acc declare create(ttlog, ddlog, yetab, helium, sica, fegroup, abartab, ebtab, wratetab, massfractab)

  private ntemp, nden, nye, npts

contains

  subroutine init_nse()

    integer :: irho, it9, iye
    integer :: j, k

    integer :: un

    ! begin initialization

    allocate(ttlog(npts))
    allocate(ddlog(npts))
    allocate(yetab(npts))
    allocate(helium(npts))
    allocate(sica(npts))
    allocate(fegroup(npts))
    allocate(abartab(npts))
    allocate(ebtab(npts))
    allocate(wratetab(npts))
    allocate(massfractab(nspec, npts))

    ! set table parameters

    ! read in table
    open(newunit=un, file="nse19.tbl")

5   format(2f12.5, 1pe12.5, 6e12.5, 19e12.5)

    do irho = 1, nden
       do it9 = 1, ntemp
          do iye = 1, nye
             j = (irho-1)*ntemp*nye + (it9-1)*nye + iye
             read (un, 5) ttlog(j), ddlog(j), yetab(j), helium(j), sica(j), &
                          fegroup(j), abartab(j), ebtab(j), wratetab(j), &
                          (massfractab(k,j), k=1, nspec)

          end do
       end do
    end do

    !$acc update device(ttlog, ddlog, yetab)
    !$acc update device(helium, sica, fegroup)
    !$acc update device(abartab, ebtab, wratetab, massfractab)

  end subroutine init_nse

  subroutine nse_interp(T, rho, ye, abar, dq, dyedt, X)

    implicit none

    real(rt), intent(in) :: T, rho, ye
    real(rt), intent(inout) :: abar, dq, dyedt
    real(rt), intent(inout) :: X(nspec)

    integer :: n

    real(rt) :: tlog, rholog, yet
    integer :: it1, it2
    integer :: ir1, ir2
    integer :: ic1, ic2

    integer :: it1r1c1
    integer :: it1r1c2
    integer :: it1r2c1
    integer :: it1r2c2
    integer :: it2r1c1
    integer :: it2r1c2
    integer :: it2r2c1
    integer :: it2r2c2

    real(rt) :: t0, r0, x0
    real(rt) :: td, rd, xd
    real(rt) :: omtd, omrd, omxd

    !$gpu

    tlog = log10(T)
    rholog = log10(rho)
    yet = ye

    if (tlog < 9.0d0) tlog = 9.0d0
    if (tlog > 10.4d0) tlog = 10.4d0

    it1 = int((tlog -9.0d0)*50.0d0 - 1.d-6)
    it1 = it1 + 1
    it2 = it1 + 1

    if (rholog < 7.0d0) rholog = 7.0d0
    if (rholog > 10.0d0) rholog = 10.0d0

    ir1 = int((rholog -7.0d0)*10.0d0 - 1.d-6)
    ir1 = ir1 + 1
    ir2 = ir1+1

    if (yet < 0.40d0) yet = 0.40d0
    if (yet > 0.50d0) yet = 0.50d0

    ic1 = int((0.50d0-yet)/0.005d0 - 1.0d-6)
    ic1 = ic1 + 1
    ic2 = ic1 + 1

    ! find the eight interpolation points in the 1D arrays

    it1r1c1 = (ir1-1)*71*21 + (it1-1)*21 + ic1
    it1r1c2 = (ir1-1)*71*21 + (it1-1)*21 + ic2
    it1r2c1 = (ir2-1)*71*21 + (it1-1)*21 + ic1
    it1r2c2 = (ir2-1)*71*21 + (it1-1)*21 + ic2
    it2r1c1 = (ir1-1)*71*21 + (it2-1)*21 + ic1
    it2r1c2 = (ir1-1)*71*21 + (it2-1)*21 + ic2
    it2r2c1 = (ir2-1)*71*21 + (it2-1)*21 + ic1
    it2r2c2 = (ir2-1)*71*21 + (it2-1)*21 + ic2

    t0 = 9.0d0 + real(it1-1)*0.02d0
    r0 = 7.0d0 + real(ir1-1)*0.10d0
    x0 = 0.50d0 - real(ic1-1)*0.005d0

    td = (tlog - t0)/0.02d0
    rd = (rholog - r0)/0.10d0
    xd = (x0-yet)/0.005d0
    xd = max(0.0d0,xd)
    omtd = 1.0d0 - td
    omrd = 1.0d0 - rd
    omxd = 1.0d0 - xd

    abar = abartab(it1r1c1)*omtd*omrd*omxd &
         + abartab(it1r1c2)*omtd*omrd*xd &
         + abartab(it1r2c1)*omtd*rd*omxd &
         + abartab(it1r2c2)*omtd*rd*xd &
         + abartab(it2r1c1)*td*omrd*omxd &
         + abartab(it2r1c2)*td*omrd*xd &
         + abartab(it2r2c1)*td*rd*omxd &
         + abartab(it2r2c2)*td*rd*xd

    dq   = ebtab(it1r1c1)*omtd*omrd*omxd &
         + ebtab(it1r1c2)*omtd*omrd*xd &
         + ebtab(it1r2c1)*omtd*rd*omxd &
         + ebtab(it1r2c2)*omtd*rd*xd &
         + ebtab(it2r1c1)*td*omrd*omxd &
         + ebtab(it2r1c2)*td*omrd*xd &
         + ebtab(it2r2c1)*td*rd*omxd &
         + ebtab(it2r2c2)*td*rd*xd

    dyedt = wratetab(it1r1c1)*omtd*omrd*omxd &
          + wratetab(it1r1c2)*omtd*omrd*xd &
          + wratetab(it1r2c1)*omtd*rd*omxd &
          + wratetab(it1r2c2)*omtd*rd*xd &
          + wratetab(it2r1c1)*td*omrd*omxd &
          + wratetab(it2r1c2)*td*omrd*xd &
          + wratetab(it2r2c1)*td*rd*omxd &
          + wratetab(it2r2c2)*td*rd*xd

    ! this is actually the sum of all e- capture and e+ decay, so if
    ! e- capture dominates, this quantity is positive, but Ye should
    ! decrease, so we swap the sign here.
    dyedt = -dyedt

    do n = 1, nspec
       X(n) = massfractab(n, it1r1c1)*omtd*omrd*omxd &
            + massfractab(n, it1r1c2)*omtd*omrd*xd &
            + massfractab(n, it1r2c1)*omtd*rd*omxd &
            + massfractab(n, it1r2c2)*omtd*rd*xd &
            + massfractab(n, it2r1c1)*td*omrd*omxd &
            + massfractab(n, it2r1c2)*td*omrd*xd &
            + massfractab(n, it2r2c1)*td*rd*omxd &
            + massfractab(n, it2r2c2)*td*rd*xd
    end do

  end subroutine nse_interp

end module nse_module

