# `xrb_ignition`

`xrb_ignition` finds the ignition condition of the X-ray burst
based on https://github.com/andrewcumming/settle
See DOI: 10.1086/317191 for details.

Usage:
./main.ex <X> <Z> <Fb> <mdot> <COMPRESS>

X: Hydrogen massfraction
Z: CNO massfraction
Fb: Base Flux in MeV/nucleon/mdot
mdot: Local accretion rate in Eddington unit
COMPRESS: To include compressional heating or not.
