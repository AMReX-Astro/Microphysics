# C-burning rate module generator

import pynucastro as pyna

def create_network():

    srates = ["c12(c12,a)ne20",
             "c12(c12,n)mg23",
             "c12(c12,p)na23",
             "c12(a,g)o16"]

    rl = pyna.ReacLibLibrary()
    rates = rl.get_rate_by_name(srates)

    tl = pyna.TabularLibrary()
    rates.append(tl.get_rate_by_name("n(,)p"))

    lib = pyna.Library(rates=rates)

    c_net = pyna.AmrexAstroCxxNetwork(libraries=[lib])
    return c_net


if __name__ == "__main__":

    c_net = create_network()
    c_net.write_network()

    c_net.summary()
