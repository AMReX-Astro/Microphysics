import pynucastro as pyna
from pynucastro import Nucleus

nuc_list = [Nucleus("n"), Nucleus("p"), Nucleus("d"),
            Nucleus("he3"), Nucleus("he4"), Nucleus("c12"), Nucleus("o16"),
            Nucleus("ne20"), Nucleus("ne21"), Nucleus("ne22"), Nucleus("n13"), Nucleus("n14"), Nucleus("f18"),
            Nucleus("na23"), Nucleus("mg24"), Nucleus("si28"),
            Nucleus("s32"), Nucleus("ar36"), Nucleus("ca40"),
            Nucleus("al27"), Nucleus("p31"), Nucleus("cl35"), Nucleus("k39")]

#nuc_list += pyna.get_nuclei_in_range(20, 20, 45, 48)
nuc_list += pyna.get_nuclei_in_range(21, 21, 43, 46)
nuc_list += pyna.get_nuclei_in_range(22, 22, 44, 49)
nuc_list += pyna.get_nuclei_in_range(23, 23, 47, 54)
nuc_list += pyna.get_nuclei_in_range(24, 24, 48, 56)
nuc_list += pyna.get_nuclei_in_range(25, 25, 51, 58)
nuc_list += pyna.get_nuclei_in_range(26, 26, 52, 60)
nuc_list += pyna.get_nuclei_in_range(27, 27, 54, 61) # 64
nuc_list += pyna.get_nuclei_in_range(28, 28, 56, 62) # 65
nuc_list.append(Nucleus("cu59"))
nuc_list.append(Nucleus("zn60"))

tl = pyna.TabularLibrary()
rl = pyna.ReacLibLibrary()
tlib = tl.linking_nuclei(nuc_list)
rlib = rl.linking_nuclei(nuc_list)

rc = pyna.RateCollection(libraries=[rlib, tlib])

dupes = rc.find_duplicate_links()

pp = dupes.pop()
from pynucastro.rates import ReacLibRate
rates_to_remove = []
for d in dupes:
    rates_to_remove += [r for r in d if isinstance(r, ReacLibRate)]
    
rc.remove_rates(rates_to_remove)

import numpy as np

Ts = np.logspace(np.log10(3.e9), np.log10(2.e10), 51)
rhos = np.logspace(7, 10, 31)
yes = np.linspace(0.43, 0.5, 15)

for T in reversed(Ts):
    initial_guess = (-3.5, -15)
    print(f"working on {T=}")
    for rho in reversed(rhos):
        for n, ye in enumerate(reversed(yes)):
            try:
                comp, sol = rc.get_comp_nse(rho, T, ye, use_coulomb_corr=True,
                                            init_guess=initial_guess, return_sol=True)
            except ValueError:
                initial_guess = (-3.5, -15)
                comp, sol = rc.get_comp_nse(rho, T, ye, use_coulomb_corr=True,
                                            init_guess=initial_guess, return_sol=True)
                
            initial_guess = sol

