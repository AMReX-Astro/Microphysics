import numpy as np
import argparse
from StarKiller.initialization import starkiller_initialize
from StarKiller.interfaces import BurnType
from StarKiller.network import Network

starkiller_initialize("probin")

burn_state = BurnType()
network = Network()

# Note that variables in the burn type are all lowercase.
# Use of the wrong case will cause a silent failure where
# the respective variable is not set.
burn_state.state.rho = 1.0e9
burn_state.state.t = 1.0e9
burn_state.state.xn = np.zeros(network.nspec)

# X(H1) = X(He4) = 0.49 with 0.02 mass fraction spread evenly throughout the network
burn_state.state.xn = 0.02/(network.nspec - 2)
burn_state.state.xn[network.species_map["h1"]] = 0.49
burn_state.state.xn[network.species_map["he4"]] = 0.49

print("Before RHS call:")
print(burn_state.state)

network.rhs(burn_state)
network.jacobian(burn_state)

print("After RHS call:")
print(burn_state.state)

def nparray_to_string(x):
    y = ["{}".format(xi) for xi in x]
    z = "  ".join(y)
    return z

f = open("rhs.txt", "w")
f.write("rp process network RHS evaluation\n")
f.write("\n density:\n {}\n".format(burn_state.state.rho))
f.write("\n temperature:\n {}\n".format(burn_state.state.t))
f.write("\n species:\n {}\n".format(nparray_to_string(network.short_species_names)))
f.write("\n X:\n {}\n".format(nparray_to_string(burn_state.state.xn)))
f.write("\n Species Ydot:\n {}\n".format(nparray_to_string(network.get_species_molar_rhs(burn_state))))
f.close()
