# C-burning with A=23 URCA rate module generator

import pynucastro as pyna

rl = pyna.ReacLibLibrary()
rl_rates = rl.get_rate_by_name(["c12(c12,a)ne20",
                                "c12(c12,n)mg23",
                                "c12(c12,p)na23",
                                "c12(a,g)o16"])

tl = pyna.TabularLibrary()
tl_rates = tl.get_rate_by_name(["na23(,)ne23",
                                "ne23(,)na23",
                                "n(,)p",
                                "p(,)n"])

urca_net = pyna.AmrexAstroCxxNetwork(rates=rl_rates+tl_rates)
urca_net.write_network()
urca_net.summary()
