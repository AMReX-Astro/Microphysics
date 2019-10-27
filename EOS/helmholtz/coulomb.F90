module coulomb_module

  use helmholtz_constants_module

  implicit none

  private

  ! Constants used for the Coulomb corrections
  double precision, parameter :: a1 = -0.898004d0
  double precision, parameter :: b1 =  0.96786d0
  double precision, parameter :: c1 =  0.220703d0
  double precision, parameter :: d1 = -0.86097d0
  double precision, parameter :: e1 =  2.5269d0
  double precision, parameter :: a2 =  0.29561d0
  double precision, parameter :: b2 =  1.9885d0
  double precision, parameter :: c2 =  0.288675d0
  double precision, parameter :: onethird = 1.0d0/3.0d0
  double precision, parameter :: esqu = qe * qe

  public :: apply_coulomb_corrections

contains

  subroutine apply_coulomb_corrections(den, temp, kt, ktinv, abar, zbar, ytot1, xni, &
                                       dpiondd, dpiondt, dxnidd, dxnida, &
                                       prad, pion, pele, erad, eion, eele, &
                                       ecoul, decouldd, decouldt, decoulda, decouldz, &
                                       pcoul, dpcouldd, dpcouldt, dpcoulda, dpcouldz, &
                                       scoul, dscouldd, dscouldt, dscoulda, dscouldz)

    use amrex_constants_module, only: ZERO

    implicit none

    double precision, intent(in   ) :: den, temp, kt, ktinv, abar, zbar, ytot1, xni
    double precision, intent(in   ) :: dpiondd, dpiondt, dxnidd, dxnida
    double precision, intent(in   ) :: prad, pion, pele
    double precision, intent(in   ) :: erad, eion, eele
    double precision, intent(inout) :: ecoul, decouldd, decouldt, decoulda, decouldz
    double precision, intent(inout) :: pcoul, dpcouldd, dpcouldt, dpcoulda, dpcouldz
    double precision, intent(inout) :: scoul, dscouldd, dscouldt, dscoulda, dscouldz

    double precision :: dsdd, dsda, lami, inv_lami, lamida, lamidd
    double precision :: plasg, plasgdd, plasgdt, plasgda, plasgdz

    double precision :: s, x, y, z
    double precision :: p_temp, e_temp

    pcoul    = ZERO
    dpcouldd = ZERO
    dpcouldt = ZERO
    dpcoulda = ZERO
    dpcouldz = ZERO
    ecoul    = ZERO
    decouldd = ZERO
    decouldt = ZERO
    decoulda = ZERO
    decouldz = ZERO
    scoul    = ZERO
    dscouldd = ZERO
    dscouldt = ZERO
    dscoulda = ZERO
    dscouldz = ZERO

    !..uniform background corrections only
    !..from yakovlev & shalybkov 1989
    !..lami is the average ion seperation
    !..plasg is the plasma coupling parameter
    z        = forth * pi
    s        = z * xni
    dsdd     = z * dxnidd
    dsda     = z * dxnida

    lami     = 1.0d0/s**onethird
    inv_lami = 1.0d0/lami
    z        = -onethird * lami
    lamidd   = z * dsdd/s
    lamida   = z * dsda/s

    plasg    = zbar*zbar*esqu*ktinv*inv_lami
    z        = -plasg * inv_lami
    plasgdd  = z * lamidd
    plasgda  = z * lamida
    plasgdt  = -plasg*ktinv * kerg
    plasgdz  = 2.0d0 * plasg/zbar

    !...yakovlev & shalybkov 1989 equations 82, 85, 86, 87
    if (plasg .ge. 1.0D0) then
       x        = plasg**(0.25d0)
       y        = avo_eos * ytot1 * kerg
       ecoul    = y * temp * (a1*plasg + b1*x + c1/x + d1)
       pcoul    = onethird * den * ecoul
       scoul    = -y * (3.0d0*b1*x - 5.0d0*c1/x &
                  + d1 * (log(plasg) - 1.0d0) - e1)

       y        = avo_eos*ytot1*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
       decouldd = y * plasgdd
       decouldt = y * plasgdt + ecoul/temp
       decoulda = y * plasgda - ecoul/abar
       decouldz = y * plasgdz

       y        = onethird * den
       dpcouldd = onethird * ecoul + y*decouldd
       dpcouldt = y * decouldt
       dpcoulda = y * decoulda
       dpcouldz = y * decouldz

       y        = -avo_eos*kerg/(abar*plasg)* &
                  (0.75d0*b1*x+1.25d0*c1/x+d1)
       dscouldd = y * plasgdd
       dscouldt = y * plasgdt
       dscoulda = y * plasgda - scoul/abar
       dscouldz = y * plasgdz

       !...yakovlev & shalybkov 1989 equations 102, 103, 104
    else if (plasg .lt. 1.0D0) then
       x        = plasg*sqrt(plasg)
       y        = plasg**b2
       z        = c2 * x - onethird * a2 * y
       pcoul    = -pion * z
       ecoul    = 3.0d0 * pcoul/den
       scoul    = -avo_eos/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

       s        = 1.5d0*c2*x/plasg - onethird*a2*b2*y/plasg
       dpcouldd = -dpiondd*z - pion*s*plasgdd
       dpcouldt = -dpiondt*z - pion*s*plasgdt
#ifdef EXTRA_THERMO
       dpcoulda = -dpionda*z - pion*s*plasgda
       dpcouldz = -dpiondz*z - pion*s*plasgdz
#endif

       s        = 3.0d0/den
       decouldd = s * dpcouldd - ecoul/den
       decouldt = s * dpcouldt
       decoulda = s * dpcoulda
       decouldz = s * dpcouldz

       s        = -avo_eos*kerg/(abar*plasg)* &
                  (1.5d0*c2*x-a2*(b2-1.0d0)*y)
       dscouldd = s * plasgdd
       dscouldt = s * plasgdt
       dscoulda = s * plasgda - scoul/abar
       dscouldz = s * plasgdz
    end if

    ! Disable Coulomb corrections if they cause
    ! the energy or pressure to go negative.

    p_temp = prad + pion + pele + pcoul
    e_temp = erad + eion + eele + ecoul

    if (p_temp .le. ZERO .or. e_temp .le. ZERO) then

       pcoul    = 0.0d0
       dpcouldd = 0.0d0
       dpcouldt = 0.0d0
       dpcoulda = 0.0d0
       dpcouldz = 0.0d0
       ecoul    = 0.0d0
       decouldd = 0.0d0
       decouldt = 0.0d0
       decoulda = 0.0d0
       decouldz = 0.0d0
       scoul    = 0.0d0
       dscouldd = 0.0d0
       dscouldt = 0.0d0
       dscoulda = 0.0d0
       dscouldz = 0.0d0

    end if
 
  end subroutine apply_coulomb_corrections

end module coulomb_module
