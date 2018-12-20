# C-burning rate module generator

from pynucastro.networks import StarKillerNetwork

files = ["c12-ag-o16-nac2"]

c_net = StarKillerNetwork(files)
c_net.write_network(use_cse=False)







