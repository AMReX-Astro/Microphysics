import pynucastro as pyna

import he_burn_core


DO_DERIVED_RATES = True


def doit():

    lib = he_burn_core.get_core_library(include_n14_approx=True,
                                        include_zn=False,
                                        include_iron_peak=True,
                                        include_low_ye=False,
                                        do_detailed_balance=DO_DERIVED_RATES)

    net = pyna.AmrexAstroCxxNetwork(libraries=[lib],
                                    symmetric_screening=False)

    # now we approximate some (alpha, p)(p, gamma) links

    net.make_ap_pg_approx(intermediate_nuclei=["cl35", "k39", "sc43", "v47"])
    net.remove_nuclei(["cl35", "k39", "sc43", "v47"])

    net.make_nn_g_approx(intermediate_nuclei=["fe53", "fe55", "ni57"])
    net.remove_nuclei(["fe53", "fe55", "ni57"])

    # make all rates with A >= 48 use NSE protons
    net.make_nse_protons(48)

    net.summary()

    # let's make a figure

    fig = net.plot(rotated=True, hide_xalpha=True,
                   size=(1800, 900),
                   node_size=600,
                   node_font_size=10)

    fig.savefig("he-burn-28anp.png", bbox_inches="tight")

    fig = net.plot(rotated=True, hide_xalpha=True,
                   size=(600, 700), Z_range=[24, 29], N_range=[-1, 4],
                   node_size=600,
                   node_font_size=10,
                   highlight_filter_function=lambda rate: isinstance(rate, pyna.rates.TabularRate) and rate.Q > 0)

    fig.savefig("he-burn-28anp-zoom.png", bbox_inches="tight")

    net.write_network()


if __name__ == "__main__":
    doit()
