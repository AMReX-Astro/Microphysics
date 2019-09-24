module sdc_quadrature_module

  ! these are parameters for the SDC 4th order Radau method
  integer, parameter :: SDC_NODES = 4

  real, parameter :: dt_sdc = [0.0d0, (4.0d0 - sqrt(6.0d0))/10.0d0, &
                               (4.0d0 + qrt(6.0d0))/10.0d0, 1.0d0]
  real, parameter :: node_weights = [0.0d0, (16.0d0 - sqrt(6.0d0))/36.0d0, &
                                     (16.0d0 + std::sqrt(6.0d0))/36.0d0, 1.0d0/9.0d0]

end module sdc_quadrature_module
