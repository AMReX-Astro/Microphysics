#Authors: Piyush Sharda & Benjamin Wibking (ANU, 2022)

# burn_cell_primordial_chem

`burn_cell_primordial_chem` integrates a primordial ISM chemistry network 
 for a single set of initial conditions.  The density, temperature, and composition 
 are set in the inputs file, as well as the maximum time to integrate.

 Upon completion, the new state is printed to the screen.

# key difference with other tests

  For primordial chemistry, state.xn is always assumed to contain
  number densities. We work with number densities and not mass fractions
  because our equations are very stiff (stiffness ratios are as high as 1e31)
  because y/ydot (natural timescale for a species abundance to vary) can be 
  very different (by factors ~ 1e30) for different species.
  However, state.rho still conatins the density in g/cm^3, and state.e 
  still contains the specific internal energy in erg/g/K.

# continuous integration

The code is built with the `primordial_chem` network and run with `inputs_primordial_chem`.
