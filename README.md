# Microphysics

*A collection of astrophysical microphysics routines with interfaces to
 the BoxLib codes*

To use this repository with BoxLib codes, set `MICROPHYSICS_DIR` to point
to the `Microphysics/` directory.

There are several core types of microphysics routines hosted here:

* `eos/`: these are the equations of state.  All of them use a Fortran
  derived type `eos_t` to pass the thermodynamic state information in
  and out.

* `networks/`: these are the reaction networks.  They serve both to
  define the composition and its properties, as well as describe the
  reactions and energy release when reactions occur.


These routines are written to be compatible with:

* Castro: http://boxlib-codes.github.io/Castro/

* Maestro: https://github.com/BoxLib-Codes/MAESTRO

