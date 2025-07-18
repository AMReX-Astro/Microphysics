import pynucastro as pyna

import he_burn_core


def create_network():

    lib = he_burn_core.get_core_library(include_n14_approx=True,
                                        include_zn=True,
                                        include_iron_peak=True,
                                        include_low_ye=True,
                                        do_detailed_balance=True)

    net = pyna.AmrexAstroCxxNetwork(libraries=[lib],
                                    symmetric_screening=False)

    # now we approximate some (alpha, p)(p, gamma) links

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47"])

    net.summary()
    return net


def doit():

    net = create_network()

    # let's make a figure

    fig = net.plot(rotated=True, hide_xalpha=True,
                   size=(1800, 900),
                   node_size=600,
                   node_font_size=10)

    fig.savefig("he-burn-33a.png", bbox_inches="tight")

    fig = net.plot(rotated=True, hide_xalpha=True,
                   size=(720, 840), Z_range=[24, 30], N_range=[-1, 4],
                   node_size=600,
                   node_font_size=10,
                   highlight_filter_function=lambda rate: isinstance(rate, pyna.rates.TabularRate))

    fig.savefig("he-burn-33a-zoom.png", bbox_inches="tight")

    net.write_network()


if __name__ == "__main__":
    doit()
