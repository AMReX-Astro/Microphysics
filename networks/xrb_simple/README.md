# xrb_simple

This is a simple network designed for capturing burning in mixed H/He
X-ray bursts.  The reaction sequence is:

3-α + 2p -> 14O (limited by the 3-α rate)

14O + α -> Ne18 (limited by 14O(α,p)17F rate)

15O + α + 6p -> 25Si (limited by 15O(α,γ)19Ne rate)

18Ne + α + 3p -> 25Si (limited by 18Ne(α,p)21Na rate)

14O + p -> 15O (limited by 14O(e+ν)14N rate)

15O + 3p -> 14O + α (limited by 15O(e+ν)15N rate)


All reactions conserve mass. Where charge is not conserved, fast weak
interactions are assumed. Weak rates are trivial, fits to the 4 strong
rates to a power law in T9=0.3−1, linear in density.
