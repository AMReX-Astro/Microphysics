import numpy as np

# needs progressbar2
import progressbar

import pynucastro as pyna
from pynucastro import Nucleus


class NSEState:
    def __init__(self, rho, T, Ye, comp, rc):
        self.rho = rho
        self.T = T
        self.Ye = Ye

        self.comp = comp

        self.ydots = rc.evaluate_ydots(self.rho, self.T, self.comp,
                                       screen_func=pyna.screening.potekhin_1998,
                                       rate_filter=lambda r: isinstance(r, pyna.rates.TabularRate))

        _, enu = rc.evaluate_energy_generation(self.rho, self.T, self.comp,
                                               screen_func=pyna.screening.potekhin_1998,
                                               return_enu=True)
        self.enu = enu

    def __str__(self):
        return f"({self.rho:12.6g}, {self.T:12.6g}, {self.Ye:6.4f}): {self.get_abar():6.4f}  {self.get_bea():6.4f}  {self.get_dyedt():12.6g}  {self.get_enu():12.6g}"

    def get_abar(self):
        return self.comp.eval_abar()

    def get_bea(self):
        return sum(q.nucbind * self.comp.X[q] for q in self.comp.X)

    def get_dyedt(self):
        return sum(q.Z * self.ydots[q] for q in self.comp.X)

    def get_dabardt(self):
        abar = self.get_abar()
        return -abar**2 * sum(self.ydots[q] for q in self.comp.X)

    def get_enu(self):
        return self.enu

    def get_aprox19_comp(self):
        aprox19_comp = [Nucleus("he3"), Nucleus("he4"), Nucleus("c12"), Nucleus("n14"),
                        Nucleus("o16"), Nucleus("ne20"), Nucleus("mg24"), Nucleus("si28"),
                        Nucleus("s32"), Nucleus("ar36"), Nucleus("ca40"), Nucleus("ti44"),
                        Nucleus("cr48"), Nucleus("fe52"), Nucleus("fe54"), Nucleus("ni56"),
                        Nucleus("n"), Nucleus("p")]
        reduced_comp = self.comp.bin_as(aprox19_comp, exclude=[Nucleus("ni56")])

        # the actual aprox19 network has 2 protons, one is the NSE proton
        # and the other is H that participates in H-burning.  We need to
        # to add H1 with zero abundance in our return

        # we also want to preserve the order that aprox19 uses
        X = []
        X.append(("H1", 0.0))
        X.append(("He3", reduced_comp.X[Nucleus("he3")]))
        X.append(("He4", reduced_comp.X[Nucleus("he4")]))
        X.append(("C12", reduced_comp.X[Nucleus("c12")]))
        X.append(("N14", reduced_comp.X[Nucleus("n14")]))
        X.append(("O16", reduced_comp.X[Nucleus("o16")]))
        X.append(("Ne20", reduced_comp.X[Nucleus("ne20")]))
        X.append(("Mg24", reduced_comp.X[Nucleus("mg24")]))
        X.append(("Si28", reduced_comp.X[Nucleus("si28")]))
        X.append(("S32", reduced_comp.X[Nucleus("s32")]))
        X.append(("Ar36", reduced_comp.X[Nucleus("ar36")]))
        X.append(("Ca40", reduced_comp.X[Nucleus("ca40")]))
        X.append(("Ti44", reduced_comp.X[Nucleus("ti44")]))
        X.append(("Cr48", reduced_comp.X[Nucleus("cr48")]))
        X.append(("Fe52", reduced_comp.X[Nucleus("fe52")]))
        X.append(("Fe54", reduced_comp.X[Nucleus("fe54")]))
        X.append(("Ni56", reduced_comp.X[Nucleus("ni56")]))
        X.append(("n", reduced_comp.X[Nucleus("n")]))
        X.append(("p", reduced_comp.X[Nucleus("p")]))

        return X

def make_nse_network():

    # list of nuclei we care about

    nuc_list = [Nucleus("n"), Nucleus("p"), Nucleus("d"),
                Nucleus("he3"), Nucleus("he4"), Nucleus("c12"), Nucleus("o16"),
                Nucleus("n13"), Nucleus("n14"), Nucleus("f18"),
                Nucleus("ne20"), Nucleus("ne21"), Nucleus("ne22"),
                Nucleus("na23"), Nucleus("mg24"), Nucleus("si28"),
                Nucleus("s32"), Nucleus("ar36"), Nucleus("ca40"), Nucleus("sc43"),
                Nucleus("al27"), Nucleus("p31"), Nucleus("cl35"), Nucleus("k39")]

    nuc_list += pyna.get_nuclei_in_range(20, 20, 45, 48)  # Ca
    nuc_list += pyna.get_nuclei_in_range(21, 21, 45, 49)  # Sc
    nuc_list += pyna.get_nuclei_in_range(22, 22, 44, 52)  # Ti
    nuc_list += pyna.get_nuclei_in_range(23, 23, 47, 54)  # V
    nuc_list += pyna.get_nuclei_in_range(24, 24, 48, 56)  # Cr
    nuc_list += pyna.get_nuclei_in_range(25, 25, 51, 58)  # Mn
    nuc_list += pyna.get_nuclei_in_range(26, 26, 52, 60)  # Fe
    nuc_list += pyna.get_nuclei_in_range(27, 27, 54, 61)  # Co
    nuc_list += pyna.get_nuclei_in_range(28, 28, 56, 65)  # Ni
    nuc_list.append(Nucleus("cu59"))
    nuc_list.append(Nucleus("zn60"))

    # create the library

    tl = pyna.TabularLibrary()
    rl = pyna.ReacLibLibrary()
    tlib = tl.linking_nuclei(nuc_list)
    rlib = rl.linking_nuclei(nuc_list)

    all_lib = rlib + tlib

    # remove dupes

    dupes = all_lib.find_duplicate_links()

    rates_to_remove = []
    for d in dupes:
        rates_to_remove += [r for r in d if isinstance(r, pyna.rates.ReacLibRate)]

    for r in rates_to_remove:
        all_lib.remove_rate(r)

    # create the rate collection

    rc = pyna.RateCollection(libraries=[all_lib])

    return rc

def generate_table():

    rc = make_nse_network()

    Ts = np.logspace(np.log10(3.e9), np.log10(2.e10), 51)
    rhos = np.logspace(7, 10, 31)
    yes = np.linspace(0.43, 0.5, 15)

    mu_p0 = -3.5
    mu_n0 = -15.0

    mu_p = np.ones((len(rhos), len(yes)), dtype=np.float64) * mu_p0
    mu_n = np.ones((len(rhos), len(yes)), dtype=np.float64) * mu_n0

    nse_states = []
    for T in reversed(Ts):
        for irho, rho in enumerate(reversed(rhos)):
            for iye, ye in enumerate(reversed(yes)):
                initial_guess = (mu_p[irho, iye], mu_n[irho, iye])
                try:
                    comp, sol = rc.get_comp_nse(rho, T, ye, use_coulomb_corr=True,
                                                init_guess=initial_guess, return_sol=True)
                except ValueError:
                    initial_guess = (-3.5, -15)
                    comp, sol = rc.get_comp_nse(rho, T, ye, use_coulomb_corr=True,
                                                init_guess=initial_guess, return_sol=True)

                mu_p[irho, iye] = sol[0]
                mu_n[irho, iye] = sol[1]

                nse_states.append(NSEState(rho, T, ye, comp, rc))
                print(nse_states[-1])

if __name__ == "__main__":
    generate_table()
