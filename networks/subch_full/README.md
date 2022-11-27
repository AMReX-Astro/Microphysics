# subch_full

This tries to recreate an aprox13 alpha-chain + including a bypass
rate for c12(a, g)o16 discussed in Shen & Bildsten (2008).

We don't approximate the rates (e.g., create an effective rate for (a,
g) and (a, p)(p, g) assuming proton equilibrium.  Therefore, we need
to explicitly include those intermediate nuclei.

Shen & Bildsten discuss the sequences:

* c14(a, g)o18(a, g)ne22 at high temperatures (T > 1 GK).  We don't consider
  this.

* n14(a, g)f18(a, p)ne21 is the one they consider important, since it
  produces protons that are then available for c12(p, g)n13(a, p)o16.

  This leaves ne21 as an endpoint, which we need to connect by
  including na22.

For the c12+c12, c12+o16, and o16+o16 rates, we also need to include
c12(c12,n)mg23(n, g)mg24, o16(o16, n)s31(n, g)s32,
o16(c12, n)si27(n, g)si28.  Since the neutron captures on those
intermediate nuclei are so fast, we leave those out and take the
forward rate to just be the first rate.  We do not include reverse
rates for these processes.

We could go further an include some iron-group nuclei as well to better
get the equilibrium state there.
