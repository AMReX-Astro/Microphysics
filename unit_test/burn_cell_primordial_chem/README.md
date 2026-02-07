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

# reference and compare workflow

Use this two-step workflow to split CPU reference generation from later
comparison runs.

1. Generate/update the saved CPU reference solution:

   ```bash
   ./main1d.gnu.ex inputs_primordial_chem
   ```

   This runs the normal Richardson-enabled path and writes
   `burn_cell_final_state.txt`.

2. Compare a future run against the saved reference:

   ```bash
   ./main1d.gnu.ex inputs_primordial_chem --compare-final-state burn_cell_final_state.txt
   ```

   In compare mode, the code does not run Richardson. It performs one batch
   run on the active backend (GPU when enabled, otherwise CPU), compares each
   zone against the saved reference, and reports per-zone `PASS`/`FAIL`.
   For `T`, `e`, and each species number density, the comparison threshold is
   `2 * |saved truncation error|` from the reference file.
