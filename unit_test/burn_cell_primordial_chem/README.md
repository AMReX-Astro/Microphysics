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
  However, state.rho still contains the density in g/cm^3, and state.e 
  still contains the specific internal energy in erg/g/K.

# continuous integration

The code is built with the `primordial_chem` network and run with `inputs_primordial_chem`.

# initial composition

The supplied inputs use a primordial total nuclei ratio D/H = 2.527e-5,
from [Cooke, Pettini & Steidel (2018)](https://arxiv.org/abs/1710.11129).
Number densities are in cm^-3. The total H and He nuclei densities remain
1.000102 and 0.0775, respectively; the total D nuclei density is 2.527257754e-5.

At the initial temperature of 100 K, H+ = 1e-4 and H2 = 1e-6 are held fixed.
The six abundances H-, H2+, D+, D-, HD+ and HD are obtained by simultaneously
balancing their chemical production and destruction rates in the primordial
network. Neutral H and D supply the remaining elemental reservoirs, and the
electron density enforces charge neutrality. He+ and He++ start at the existing
1e-100 numerical floor because their ionization sources are negligible at
100 K. This is a partial steady state, not full chemical or thermal equilibrium;
it avoids artificial initial transients caused by assigning all trace species
an arbitrary 1e-40 abundance. The balance applies to these inputs at 100 K and
must be recomputed if the temperature, rates or reservoir abundances change.

Species tolerances are rtol = 1e-6 and atol = 1e-20 cm^-3, with energy
rtol = 1e-8 and atol = 1e-6. The previous species absolute tolerance of
1e-4 cm^-3 exceeded the entire initial deuterium reservoir. The stored
`reference_solution.out` and `state_over_time.txt` correspond to the supplied
inputs with the default VODE integrator.
