import pynucastro as pyna


def create_network(network_type="amrex"):

    nuc = ["n", "h1", "he4", "c12", "o16",
           "ne20", "na23", "mg23-24", "al27",
           "si27-31", "p29,31,33", "s32-35", "cl33,35,37",
           "ar36-38", "k37,39", "ca40-42", "sc41-45",
           "ti44-48", "v47-51", "cr48-54", "mn51-55",
           "fe52-58", "co55-59", "ni56-60", "cu59", "zn60"]

    net = pyna.network_helper(nuc, network_type=network_type)

    intermediate_n_nuc = ["ni57", "ni59", "co58", "fe53", "fe55",
                          "fe57", "mn52", "mn54", "cr49", "cr51",
                          "cr53", "v48", "v50", "ti45", "ti47", 
                          "sc42", "sc44", "ca41"]

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
