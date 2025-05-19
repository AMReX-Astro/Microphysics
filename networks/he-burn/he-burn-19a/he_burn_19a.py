import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

import he_burn_core


DO_DERIVED_RATES = True


def doit():

    subch = he_burn_core.get_core_library(include_n14_approx=True,
                                          include_zn=False,
                                          do_detailed_balance=DO_DERIVED_RATES)

    # these are the rates that we are going to allow to be optionally
    # zeroed
    r1 = subch.get_rate_by_name("c12(p,g)n13")
    r2 = subch.get_rate_by_name("n13(he4,p)o16")

    net = AmrexAstroCxxNetwork(libraries=[subch],
                               symmetric_screening=False,
                               disable_rate_params=[r1, r2])

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

    net.summary()

    net.plot(outfile="he-burn-19a.png",
             rotated=True, hide_xalpha=True,
             size=(1500, 450),
             node_size=600, node_font_size=11,
             Z_range=(1, 29))

    net.write_network()


if __name__ == "__main__":
    doit()
