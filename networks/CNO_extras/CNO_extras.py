import pynucastro as pyna


def create_network():

    net = pyna.common_networks.cno(network_type="amrex")
    net.add_inert_nucleus("fe56")

    return net


def doit():

    net = create_network()

    net.write_network()

    comp = pyna.Composition(net.get_nuclei(), init="solar")

    rho = 1.e4
    T = 1.e8

    net.plot(rho, T, comp, outfile="cno_extras.png",
             node_size=500, node_font_size="9",
             Z_range=[0, 13], N_range=[-0.5, 12.5],
             always_show_p=True)

    net.summary()

if __name__ == "__main__":
    doit()
