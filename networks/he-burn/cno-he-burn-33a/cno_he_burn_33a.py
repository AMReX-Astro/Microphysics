# a blend of CNO_extras and subch_simple

import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

import he_burn_core


DO_DERIVED_RATES = True


def doit():

    extra_reactants = ["c13", "n14", "n15",
                       "o14", "o15", "o17", "o18",
                       "f17", "f18", "f19",
                       "ne18", "ne19",
                       "ne21", "na22", "na23",
                       "mg22"]

    subch = he_burn_core.get_core_library(include_n14_approx=False,
                                          include_zn=False,
                                          extra_nuclei=extra_reactants,
                                          do_detailed_balance=DO_DERIVED_RATES)

    net = AmrexAstroCxxNetwork(libraries=[subch], symmetric_screening=False)

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43",
                                               "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

    net.summary()

    comp = pyna.Composition(net.get_nuclei())
    comp.set_all(0.1)
    comp.set_nuc("he4", 0.95)
    comp.normalize()

    rho = 1.e6
    T = 1.e9

    fig = net.plot(rotated=True, hide_xalpha=True,
                   size=(1800, 600),
                   node_size=600, node_font_size=10,
                   Z_range=(1, 29))

    net.write_network()

    fig.savefig("cno-he-burn-33a.png", bbox_inches="tight")


if __name__ == "__main__":
    doit()
