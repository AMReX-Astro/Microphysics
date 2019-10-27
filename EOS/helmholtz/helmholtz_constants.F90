module helmholtz_constants_module

  implicit none

  ! 2006 CODATA physical constants

  ! Math constants
  double precision, parameter :: pi       = 3.1415926535897932384d0

  ! Physical constants
  double precision, parameter :: h       = 6.6260689633d-27
  double precision, parameter :: qe      = 4.8032042712d-10
  double precision, parameter :: avo_eos = 6.0221417930d23
  double precision, parameter :: clight  = 2.99792458d10
  double precision, parameter :: kerg    = 1.380650424d-16
  double precision, parameter :: amu     = 1.66053878283d-24

#ifdef RADIATION
  double precision, parameter :: ssol    = 0.0d0
#else
  double precision, parameter :: ssol    = 5.67051d-5
#endif
  double precision, parameter :: asol    = 4.0d0 * ssol / clight

  ! Some other useful combinations of the constants
  double precision, parameter :: sioncon = (2.0d0 * pi * amu * kerg)/(h*h)
  double precision, parameter :: forth   = 4.0d0/3.0d0
  double precision, parameter :: kergavo = kerg * avo_eos
  double precision, parameter :: asoli3  = asol/3.0d0
  double precision, parameter :: light2  = clight * clight

end module helmholtz_constants_module
