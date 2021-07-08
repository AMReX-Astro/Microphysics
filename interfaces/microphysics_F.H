extern "C" {

  void microphysics_initialize(const amrex::Real small_temp, const amrex::Real small_dens);

  void microphysics_finalize();
}
