import pynucastro as pyna


def create_network(network_type="amrex"):

    nuc = ["n", "h1", "he4", "c12", "o16",
           "ne20", "na23", "mg23-24", "al27",
           "si27-30", "p29,31", "s30-34", "cl33,35",
           "ar34-38", "k37,39", "ca38-42", "sc41-43",
           "ti42-46", "v45-49", "cr46-52", "mn49-53",
           "fe50-56", "co53-57", "ni56-58", "cu59", "zn60"]

    net = pyna.network_helper(nuc, network_type=network_type)

    intermediate_n_nuc = ["ni57", "co54", "fe51", "fe53", "fe55",
                          "mn50", "mn52", "cr47", "cr49", "cr51",
                          "v46", "v48", "ti43", "ti45", "sc42",
                          "ca39", "ca41", "ar35", "ar37",
                          "si29", "s31", "s33"]

    net.make_nn_g_approx(intermediate_nuclei=intermediate_n_nuc)
    net.remove_nuclei(intermediate_n_nuc)

    net.make_CO_burning_approx("C")
    net.remove_nuclei(["na23"])

    return net


def doit():

    net = create_network()

    net.summary()

    fig = net.plot(rotated=True, size=(1300, 700),
                   node_size=500, node_font_size="9",
                   legend_coord=(4, 3),
                   highlight_filter_function=lambda r: isinstance(r, pyna.rates.TabularWeakRate))

    fig.savefig("massive-star-medium.png")

    net.write_network()


if __name__ == "__main__":
    doit()
