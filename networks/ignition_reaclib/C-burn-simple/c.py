# C-burning rate module generator

from pynucastro.networks import AmrexAstroCxxNetwork

files = ["c12-c12a-ne20-cf88",
         "c12-c12n-mg23-cf88",
         "c12-c12p-na23-cf88",
         "c12-ag-o16-nac2",
         "n--p-wc12"]

c_net = AmrexAstroCxxNetwork(files)
c_net.write_network()



