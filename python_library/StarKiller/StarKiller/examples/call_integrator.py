import numpy as np
import argparse
from StarKiller.initialization import starkiller_initialize
from StarKiller.interfaces import BurnType
from StarKiller.integration import Integrator
from StarKiller.network import Network

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--density', type=float, default=1.0e9,
                    help='Density in g/cm^3')
parser.add_argument('-t', '--temperature', type=float, default=1.0e8,
                    help='Temperature in K')
parser.add_argument('-x', '--mass_fractions', type=float, nargs='*',
                    help='Mass fractions of species. If not supplied, uniform mass fractions are used.')
parser.add_argument('-dt', '--time', type=float, default=1.0e-6,
                    help='Time to integrate in seconds.')
parser.add_argument('-pinit', '--probin_initialize', type=str, default='probin_aprox13',
                    help='Probin file to use for initialization.')
args = parser.parse_args()

starkiller_initialize(args.probin_initialize)

net = Network()
integrator = Integrator()

state = BurnType()

# Note that variables in the burn type are all lowercase.
# Use of the wrong case will cause a silent failure where
# the respective variable is not set.
state.state.rho = args.density
state.state.t = args.temperature

if args.mass_fractions:
    state.state.xn = np.array(args.mass_fractions)
else:
    state.state.xn = np.ones(net.nspec) / net.nspec

print("Before integrator call:")
print(state.state)

state = integrator.integrate(state, args.time)

print("After integrator call:")
print(state.state)
