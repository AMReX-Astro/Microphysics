import pynucastro as pyna


def create_network():

    net = pyna.network_helper(["h1", "he4",
                               "c12", "c13",
                               "n13", "n14", "n15",
                               "o14", "o15", "o16", "o17", "o18",
                               "f17", "f18", "f19",
                               "ne18", "ne19", "ne20",
                               "mg22", "mg24"],
                              tabular_ordering=["ffn", "langanke", "oda"],
                              inert_nuclei=["fe56"], network_type="amrex")

    return net


def doit():

    net = create_network()

    net.write_network()

    comp = pyna.Composition(net.get_nuclei())
    comp.set_solar_like()

    rho = 1.e6
    T = 1.e8

    net.plot(rho, T, comp, outfile="cno_extras.png",
             Z_range=[1, 13], N_range=[1, 13])

    net.plot(outfile="cno_extras_hide_alpha.png",
             Z_range=[1, 13], N_range=[1, 13],
             rotated=True,
             hide_xalpha=True)


if __name__ == "__main__":
    doit()
