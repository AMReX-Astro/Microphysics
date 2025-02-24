import pynucastro as pyna
from pynucastro.networks import AmrexAstroCxxNetwork

import he_burn_core


DO_DERIVED_RATES = True


def doit():

    lib = he_burn_core.get_core_library(include_n14_sequence=True,
                                        include_zn=True,
                                        include_iron_peak=True,
                                        include_low_ye=True,
                                        do_detailed_balance=DO_DERIVED_RATES)

    net = pyna.AmrexAstroCxxNetwork(libraries=[lib],
                                    symmetric_screening=False)

    # now we approximate some (alpha, p)(p, gamma) links

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47"])

    print(f"number of nuclei = {len(net.unique_nuclei)}")
    print(f"number of ReacLib rates = {len(net.reaclib_rates)}")
    print(f"number of tabular rates = {len(net.tabular_rates)}")

    # let's make a figure

    comp = pyna.Composition(net.unique_nuclei)
    comp.set_equal()

    rho = 9.e7
    T = 6.e9

    fig = net.plot(rho, T, comp,
                   rotated=True, curved_edges=True, hide_xalpha=True,
                   size=(1800, 900),
                   node_size=500, node_shape="s", node_color="#337dff",
                   node_font_size=10)

    fig.savefig("he-burn-36a.png", bbox_inches="tight")

    fig = net.plot(rotated=True, curved_edges=True, hide_xalpha=True,
                   size=(720, 840), Z_range=[24, 30], N_range=[-1, 4],
                   node_size=500, node_shape="s", node_color="#337dff",
                   node_font_size=10,
                   highlight_filter_function=lambda rate: isinstance(rate, pyna.rates.TabularRate))

    fig.savefig("he-burn-36a-zoom.png", bbox_inches="tight")

    net.write_network()


if __name__ == "__main__":
    doit()
