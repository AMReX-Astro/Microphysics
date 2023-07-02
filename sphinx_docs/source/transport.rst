**********************
Transport Coefficients
**********************

Thermal Conductivity
====================

Thermal conductivities are provided by the conductivity/
directory. At the moment, there is a single version,
stellar [1]_ that is useful
for stellar interiors.

**Important: it is assumed that the state is thermodynamically consistent
before calling the conductivity routine.** It may be necessary to do an EOS
call first, to enforce the consistency.

.. [1]
   this code comes from Frank Timmes’ website,
   https://cococubed.com/code_pages/kap.shtml
