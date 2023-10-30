# Tabular NSE

``nse_tabular`` provides support for reading in and interpolating from
a tabulation of an NSE state from a large collection of nuclei.  The
idea is that this table will be used where the thermodynamic state has
entered NSE (high T, rho) and a regular network will be used
elsewhere.

This requires compiling with

```
USE_TABULAR_NSE=TRUE
```

This will change the equation of state to work in terms of (Y_e, abar,
B/A) and those quantities will be added to the aux state in the
network.

Interpolation from the table is done with a tricubic interpolating
polynomial.


## Table contents

The table provides:

* log10(rho) : density in CGS

* log10(T) : temperature in K

* Ye : electron fraction

* Abar : mean molecular weight ($1/\bar{A} = \sum_k X_k / A_k$)

* <B/A> : average binding energy per nucleon in MeV
 ($\langle B/A \rangle = \sum_k X_k B_k / A_k$)

* dYe/dt : evolution of the electron fraction in 1/s

* dAbar/dt : evolution of Abar in 1/s

* e_nu : weak rate neutrino loss energy in erg/g/s

* X(A), X(B), ... : the reduced composition mass fractions.  They are
  assumed to be in the same order as the nuclei in the on-grid network.


## Generating the table

Table generation is managed via pynucastro.  For the current table
in Microphysics, the script `make_nse_table.py` will create the
NSE state at a number of table grid points and output the table
as a file.  There are a few things to control here:

* The nuclei used in the NSE calculation.

  Currently we use 98.  The main limits are lake of weak rates
  for some with lower A and lack of reliable spins for very
  neutron rich nuclei.

* The "on-grid" network to reduce down to.

  This requires writing a function that bins the NSE nuclei down to
  the nuclei carried on the grid.  Presently the script does this for
  `aprox19`.

* The number of grid points and the extrema

  It doesn't make sense to consider temperatures below 1.e9 K and very
  low Ye requires a lot of nuclei with low Z/A to ensure that the NSE
  solver converges well.

  Sometimes the solver can fail to converge, but may do so if a better
  initial guess for the chemical potentials is provided.  There is an
  attempt to cache the chemical potentials that worked for the last
  temperature to hopefully accelerate the convergence.

