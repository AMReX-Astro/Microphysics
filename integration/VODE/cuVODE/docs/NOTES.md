# Notes on the CUDA port

## CONDOPT revisions

The condopt array (now in rwork_t) has been shortened to 4 elements
and unused elements deleted.

The mapping from prior elements in condopt to their new indices is:

1 => 1
5 => 2
6 => 3
7 => 4

## RWORK summary

Here's what rwork contains:

- CONDOPT: Conditional or optional input/output arguments to VODE
- YH: The Nordsieck array of length neq * (maximum order + 1). We use maximum order = 5.
- WM: Array of length 2 + 2*neq**2
  - WM(1) = sqrt(UROUND)
  - WM(2) = HRL1 (old value of H*RL1 used if MITER=3)
  - WM(3:3 + neq**2 - 1) = Matrix P = I - H*RL1*J
  - WM(3 + neq**2:2 + 2*neq**2) = Saved Jacobian matrix J
- EWT: Error weight vector of length neq
- SAVF: RHS array of length neq
- ACOR: Correction array of length neq

So all told, aside from the arguments in CONDOPT, the size of the work
array scales with neq as:

L = 9*neq + 2*neq**2

The dvode subroutine copies a section of WM to SAVF when NQ is less
than MAXORD. I might need to revise that if I don't store the Jacobian
in WM. Cf. dvode.F90 line 305.

So first, I need to revise the network interface to allow evaluating
components of Jac individually.

Then, eliminate the Jacobian and P matrix storage and revise the
linear algebra routines to compute the P elements as needed.
