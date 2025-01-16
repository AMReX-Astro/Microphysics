# a blend of CNO_extras and subch_simple

import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

import he_burn_core


DO_DERIVED_RATES = False


def doit():

    extra_reactants = ["c13", "n15",
                       "o14", "o15", "o17", "o18",
                       "f17", "f19",
                       "ne18", "ne19",
                       "na22", "na23",
                       "mg22"]

    subch = he_burn_core.get_core_library(include_n14_sequence=True,
                                          include_zn=False,
                                          extra_nuclei=extra_reactants,
                                          do_detailed_balance=DO_DERIVED_RATES)

    net = AmrexAstroCxxNetwork(libraries=[subch], symmetric_screening=True)

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

    print(f"number of nuclei: {len(net.unique_nuclei)}")
    print(f"number of rates: {len(net.rates)}")

    comp = pyna.Composition(net.get_nuclei())
    comp.set_all(0.1)
    comp.set_nuc("he4", 0.95)
    comp.normalize()

    rho = 1.e6
    T = 1.e9

    net.plot(rho, T, comp, outfile="cno-he-burn-33a.png",
             rotated=True, hide_xalpha=True, curved_edges=True,
             size=(1500, 450),
             node_size=500, node_font_size=11, node_color="#337dff", node_shape="s",
             Z_range=(1, 29))

    net.write_network()


if __name__ == "__main__":
    doit()
