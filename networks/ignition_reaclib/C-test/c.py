# C-burning rate module generator

from pyreaclib.networks import BoxLibNetwork

files = ["c12-ag-o16-nac2"]

c_net = BoxLibNetwork(files, use_cse=True)
c_net.write_network()




