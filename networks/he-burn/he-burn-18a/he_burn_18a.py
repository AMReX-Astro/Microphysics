# Network for He/C burning with key rates
# to bypass the C12(a,g)O16 rate.

import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

import he_burn_core


DO_DERIVED_RATES = True


def doit():

    subch = he_burn_core.get_core_library(include_n14_sequence=False,
                                          include_zn=False,
                                          do_detailed_balance=DO_DERIVED_RATES)

    # these are the rates that we are going to allow to be optionally
    # zeroed
    r1 = subch.get_rate_by_name("c12(p,g)n13")
    r2 = subch.get_rate_by_name("n13(he4,p)o16")

    net = AmrexAstroCxxNetwork(libraries=[subch], symmetric_screening=False, disable_rate_params=[r1, r2])
    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47", "mn51", "co55"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47", "mn51", "co55"])

    # finally, the aprox nets don't include the reverse rates for
    # C12+C12, C12+O16, and O16+O16, so remove those

    print(f"number of nuclei: {len(net.unique_nuclei)}")
    print(f"number of rates: {len(net.rates)}")

    comp = pyna.Composition(net.get_nuclei())
    comp.set_all(0.1)
    comp.set_nuc("he4", 0.95)
    comp.normalize()

    rho = 1.e6
    T = 1.e9

    net.plot(rho, T, comp, outfile="he-burn-18a.png",
             rotated=True, hide_xalpha=True, curved_edges=True,
             size=(1500, 450),
             node_size=500, node_font_size=11, node_color="#337dff", node_shape="s",
             Z_range=(1, 29))

    net.write_network()


if __name__ == "__main__":
    doit()
