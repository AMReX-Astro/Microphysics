# a blend of CNO_extras and subch_simple

import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

import he_burn_core


def create_network():

    extra_reactants = ["c13", "n15",
                       "o14", "o15", "o17", "o18",
                       "f17", "f19",
                       "ne18", "ne19", "ne22",
                       "na21", "na22", "na23",
                       "al25", "al26", "p29", "p30",
                       "cl33", "cl34", "k37", "k38",
                       "sc41", "sc42", "v45", "v46",
                       "mn49", "mn50", "co54"]

    subch = he_burn_core.get_core_library(include_n14_sequence=True,
                                          remove_al27_alpha_links=False,
                                          include_zn=False,
                                          include_iron_peak=True,
                                          extra_nuclei=extra_reactants,
                                          do_detailed_balance=True)

    net = AmrexAstroCxxNetwork(libraries=[subch])

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47"])

    return net

def doit():

    net = create_network()

    net.summary()

    fig = net.plot(rotated=True, hide_xalpha=True, curved_edges=False,
                   size=(1500, 450),
                   node_size=500, node_font_size=10,
                   Z_range=(1, 29),
                   highlight_filter_function=lambda r : isinstance(r, pyna.rates.TabularRate))

    net.write_network()

    fig.savefig("cno-he-burn-big.png", bbox_inches="tight")

if __name__ == "__main__":
    doit()
