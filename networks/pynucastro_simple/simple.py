import pynucastro as pyna
rl = pyna.ReacLibLibrary()
lib = rl.linking_nuclei(["he4", "c12", "o16", "ne20"])
net = pyna.StarKillerCxxNetwork(libraries=[lib])
net.write_network()
