# ``screen_derivs_check``

This is a simple unit test for the screening routines full set of
derivatives: simultaneous T, M₁ = Σᵢ (Zᵢ Yᵢ), and M₂ = Σᵢ (Zᵢ² Yᵢ).
We do a simple centered-difference and compare to the function call /
autodiff for a few different states.

It should be repeat with ``SCREEN_METHOD=screen5``, ``chugunov2007``,
``chugunov2009``, and ``chabrier1998``.

The Chugunov 2009 and Chabrier 1998 checks use stronger coupling and a
0.5% derivative tolerance because ``fast_atan``'s existing autodiff rule uses
the exact atan derivative while its value uses an approximation. Other
methods use a 2e-5 relative derivative tolerance, with a roundoff floor.
