import numpy as np

import pynucastro as pyna
from pynucastro import Nucleus


def get_aprox19_comp(comp):
    aprox19_comp = [Nucleus("he3"), Nucleus("he4"), Nucleus("c12"), Nucleus("n14"),
                    Nucleus("o16"), Nucleus("ne20"), Nucleus("mg24"), Nucleus("si28"),
                    Nucleus("s32"), Nucleus("ar36"), Nucleus("ca40"), Nucleus("ti44"),
                    Nucleus("cr48"), Nucleus("fe52"), Nucleus("fe54"), Nucleus("ni56"),
                    Nucleus("n"), Nucleus("p")]
    reduced_comp = comp.bin_as(aprox19_comp, exclude=[Nucleus("ni56")])

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

    rc = pyna.NSENetwork(libraries=[all_lib])

    return rc

def output_header(Ts, rhos, yes):

    with open("nse_table_size.H", "w") as nse_h:

        nse_h.write("#ifndef NSE_TABLE_SIZE_H\n")
        nse_h.write("#define NSE_TABLE_SIZE_H\n\n")

        nse_h.write("namespace nse_table_size {\n\n")

        nse_h.write(f"    constexpr int ntemp = {len(Ts)};\n")
        nse_h.write(f"    constexpr int nden = {len(rhos)};\n")
        nse_h.write(f"    constexpr int nye = {len(yes)};\n\n")

        nse_h.write(f"    constexpr Real logT_min = {np.log10(Ts.min())};\n")
        nse_h.write(f"    constexpr Real logT_max = {np.log10(Ts.max())};\n")
        nse_h.write(f"    constexpr Real dlogT = {(np.log10(Ts.max()) - np.log10(Ts.min())) / (len(Ts) - 1)}\n\n")

        nse_h.write(f"    constexpr Real logrho_min = {np.log10(rhos.min())};\n")
        nse_h.write(f"    constexpr Real logrhon_max = {np.log10(rhos.max())};\n")
        nse_h.write(f"    constexpr Real dlogrho = {(np.log10(rhos.max()) - np.log10(rhos.min())) / (len(rhos) - 1)}\n\n")

        nse_h.write(f"    constexpr Real ye_min = {yes.min()};\n")
        nse_h.write(f"    constexpr Real ye_max = {yes.max()};\n")
        nse_h.write(f"    constexpr Real dye = {(yes.max() - yes.min()) / (len(yes) - 1)}\n\n")

        nse_h.write("}\n")
        nse_h.write("#endif\n")


def generate_table():

    nse_net = make_nse_network()

    Ts = np.logspace(9.4, 10.4, 51)
    rhos = np.logspace(7, 10, 31)
    yes = np.linspace(0.43, 0.5, 15)

    output_header(Ts, rhos, yes)

    mu_p0 = -3.5
    mu_n0 = -15.0

    mu_p = np.ones((len(rhos), len(yes)), dtype=np.float64) * mu_p0
    mu_n = np.ones((len(rhos), len(yes)), dtype=np.float64) * mu_n0

    nse_net.generate_table(rho_values=rhos,
                           T_values=Ts,
                           Ye_values=yes,
                           comp_reduction_func=get_aprox19_comp,
                           verbose=True)

if __name__ == "__main__":
    generate_table()
