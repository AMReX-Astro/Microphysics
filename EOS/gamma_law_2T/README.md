# gamma_law_2T

This is an equation of state for a 2-temperature fluid.  We assume
that there is a composition that has neutrals and ions, with
effectively the same mass.  We call these collectively the "heavies".
The electrons are assumed to have zero mass but still contribute to
the thermodynamics.

We assume that there is a single gamma for all components of the EOS.

At the moment, we assume only a single neutral (this has Z = 0) and a
single ion (this has Z > 0).  But this can be generalized as needed.

This also requires 2 pieces of auxiliary data:

   * the specific internal energy of all the "heavies"
   * the specific internal energy of the electrons.

It is suggested that you use the `general_null` network with the
`gammalaw_2T.net` inputs file.

There is really no single temperature, so inputs like `eos_input_re`
that want to return a unique temperature will instead get an effective
average temperature, such that calling the EOS with this average
temperature will give the same input e.

