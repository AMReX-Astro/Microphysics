An explicit integrator using the quasi-steady-state approximation as discussed in

M W Guidry and J A Harris, Computational Science & Discovery 6 (2013) 015002,
"Explicit integration of extremely stiff reaction networks: quasi-steady-state methods"

and originally proposed in

David R. Mott, Elaine S. Oran, Bram van Leer, JCP 164, 407 (2000),
"A Quasi-Steady-State Solver for the Stiff Ordinary Differential Equations of Reaction Kinetics"

The timestepping is determined by not letting any integration
quantity change by more than a given factor in a timestep.
