import numpy as np
import argparse
from StarKiller.initialization import starkiller_initialize
from StarKiller.interfaces import EosType
from StarKiller.eos import Eos

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--density', type=float, default=1.0e9,
                    help='Density in g/cm^3')
parser.add_argument('-t', '--temperature', type=float, default=1.0e8,
                    help='Temperature in K')
parser.add_argument('-x', '--mass_fractions', type=float, nargs='*',
                    help='Mass fractions of species. Overrides --abar and --zbar if supplied.')
parser.add_argument('-a', '--abar', type=float, default=12.0,
                    help='Average atomic mass A')
parser.add_argument('-z', '--zbar', type=float, default=6.0,
                    help='Average atomic number Z')
parser.add_argument('-pinit', '--probin_initialize', type=str, default='probin_aprox13',
                    help='Probin file to use for initialization.')
args = parser.parse_args()

starkiller_initialize(args.probin_initialize)

eos_state = EosType()

# Note that variables in the eos type are all lowercase.
# Use of the wrong case will cause a silent failure where
# the respective variable is not set.
eos_state.state.rho = args.density
eos_state.state.t = args.temperature

use_raw_inputs = False
if args.mass_fractions:
    eos_state.state.xn = np.array(args.mass_fractions)
    use_raw_inputs = False
else:
    eos_state.state.abar = args.abar
    eos_state.state.zbar = args.zbar
    eos_state.state.y_e = args.zbar/args.abar
    eos_state.state.mu_e = 1.0/eos_state.state.y_e
    use_raw_inputs = True

print("Before EOS call:")
print(eos_state.state)

# use_raw_inputs = True will tell the EOS interface not to
# use the mass fractions to set abar, zbar before calling
# the actual EOS subroutine.
eos = Eos()
eos.evaluate(EosType.eos_input_rt, eos_state, use_raw_inputs=use_raw_inputs)

print("After EOS call:")
print(eos_state.state)
