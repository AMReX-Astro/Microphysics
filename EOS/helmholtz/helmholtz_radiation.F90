module helmholtz_radiation_module

  use helmholtz_constants_module, only: asoli3

  implicit none

contains

  subroutine apply_radiation(deni, temp, tempi, &
                             prad, dpraddd, dpraddt, dpradda, dpraddz, &
                             erad, deraddd, deraddt, deradda, deraddz, &
                             srad, dsraddd, dsraddt, dsradda, dsraddz)

    implicit none

    double precision, intent(in   ) :: deni, temp, tempi
    double precision, intent(inout) :: prad, dpraddd, dpraddt, dpradda, dpraddz
    double precision, intent(inout) :: erad, deraddd, deraddt, deradda, deraddz
    double precision, intent(inout) :: srad, dsraddd, dsraddt, dsradda, dsraddz

    prad    = asoli3 * temp * temp * temp * temp
    dpraddd = 0.0d0
    dpraddt = 4.0d0 * prad*tempi
#ifdef EXTRA_THERMO
    dpradda = 0.0d0
    dpraddz = 0.0d0
#endif

    erad    = 3.0d0 * prad*deni
    deraddd = -erad*deni
    deraddt = 3.0d0 * dpraddt*deni
#ifdef EXTRA_THERMO
    deradda = 0.0d0
    deraddz = 0.0d0
#endif

    srad    = (prad*deni + erad)*tempi
    dsraddd = (dpraddd*deni - prad*deni*deni + deraddd)*tempi
    dsraddt = (dpraddt*deni + deraddt - srad)*tempi
#ifdef EXTRA_THERMO
    dsradda = 0.0d0
    dsraddz = 0.0d0
#endif

  end subroutine apply_radiation

end module helmholtz_radiation_module
