*******
Preface
*******

Welcome to the AMReX-Astro Microphysics!

In this User’s Guide we describe the microphysics modules designed to
enable simulations of stellar explosions.

The original design was to support the AMReX codes CASTRO and
MAESTRO. These all have a consistent interface and are designed to
provide the users of those codes an easy experience in moving from the
barebones microphysics modules provided in those codes. For the
purposes of this user’s guide, the microphysical components we
currently deal with are the equation of state (EOS) and the nuclear
burning network.

Microphysics is not a stand-alone code. It is intended to be used in
conjunction with a simulation code. At the moment, the interfaces and
build stubs are compatible with the AMReX codes. In many cases we
will provide test modules that demonstrate a minimal working example
for how to run the modules (leveraging the AMReX build system). The
goal is to make the support more general, and extend to other codes
in the future.

A number of the routines contained here we authored by other people.
We bundle them here with permission, usually changing the interfaces
to be compatible with our standardized interface. We in particular
thank Frank Timmes for numerous reaction networks and his equation
of state routines.
