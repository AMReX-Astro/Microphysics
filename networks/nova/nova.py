import pynucastro as pyna


def create_network():

    all_nuclei = ["p", "h2", "he3", "he4", "be7", "b8",
                  "c12", "c13", "n13", "n14", "n15",
                  "o14", "o15", "o16", "o17", "f17", "f18"]

    net = pyna.network_helper(all_nuclei,
                              tabular_ordering=["ffn", "langanke", "oda"],
                              network_type="amrex")

    return net


if __name__ == "__main__":

    net = create_network()

    net.write_network()

    comp = pyna.Composition(net.get_nuclei())
    comp.set_solar_like()

    rho = 1.e3
    T = 1.e8

    edge_labels = {(pyna.Nucleus("he4"), pyna.Nucleus("c12")):
                   r"$\alpha(\alpha\alpha,\gamma){}^{12}\mathrm{C}$"}

    net.plot(rotated=False, outfile="nova.png",
             N_range=(-1, 10), Z_range=(0, 10),
             hide_xalpha=True, hide_xp=True,
             edge_labels=edge_labels,
             node_size=300, node_font_size=10)

    net.summary()
