# rprox

This is an approximate network for rp-process burning.  It is
essentially the network described in Appendix C of Wallace & Woosley,
ApJS 45, 389 (1981) but with updated reaction rates from ReacLib.

The default setup is for burning in a neutron star atmosphere.  If you
want a setup at lower density, you need to change the Lweak paramter
in actual_rhs.f90.

This network was used in:

 * Malone et al. 2014
   https://ui.adsabs.harvard.edu/abs/2014ApJ...788..115M/abstract

 * Zingale et al. 2015
   https://ui.adsabs.harvard.edu/abs/2015ApJ...807...60Z/abstract

