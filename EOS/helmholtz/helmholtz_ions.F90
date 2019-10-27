module helmholtz_ions_module

  use helmholtz_constants_module, only: avo_eos, kerg, kergavo, sioncon

  implicit none

contains

  subroutine apply_ions(den, deni, temp, tempi, kt, abar, ytot1, &
                        xni, dxnidd, dxnida, &
                        pion, dpiondd, dpiondt, dpionda, dpiondz, &
                        eion, deiondd, deiondt, deionda, deiondz, &
                        sion, dsiondd, dsiondt, dsionda, dsiondz)

    implicit none

    double precision, intent(in   ) :: den, deni, temp, tempi, kt, abar, ytot1
    double precision, intent(inout) :: xni, dxnidd, dxnida
    double precision, intent(inout) :: pion, dpiondd, dpiondt, dpionda, dpiondz
    double precision, intent(inout) :: eion, deiondd, deiondt, deionda, deiondz
    double precision, intent(inout) :: sion, dsiondd, dsiondt, dsionda, dsiondz

    double precision :: s, x, y, z

    xni     = avo_eos * ytot1 * den
    dxnidd  = avo_eos * ytot1
    dxnida  = -xni * ytot1

    pion    = xni * kt
    dpiondd = dxnidd * kt
    dpiondt = xni * kerg
#ifdef EXTRA_THERMO
    dpionda = dxnida * kt
    dpiondz = 0.0d0
#endif

    eion    = 1.5d0 * pion*deni
    deiondd = (1.5d0 * dpiondd - eion)*deni
    deiondt = 1.5d0 * dpiondt*deni
#ifdef EXTRA_THERMO
    deionda = 1.5d0 * dpionda*deni
    deiondz = 0.0d0
#endif

    x       = abar*abar*sqrt(abar) * deni/avo_eos
    s       = sioncon * temp
    z       = x * s * sqrt(s)
    y       = log(z)
    sion    = (pion*deni + eion)*tempi + kergavo * ytot1 * y
    dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi &
              - kergavo * deni * ytot1
    dsiondt = (dpiondt*deni + deiondt)*tempi -  &
              (pion*deni + eion) * tempi*tempi  &
              + 1.5d0 * kergavo * tempi*ytot1
    x       = avo_eos*kerg/abar
#ifdef EXTRA_THERMO
    dsionda = (dpionda*deni + deionda)*tempi  &
              + kergavo*ytot1*ytot1* (2.5d0 - y)
    dsiondz = 0.0d0
#endif

  end subroutine apply_ions

end module helmholtz_ions_module
