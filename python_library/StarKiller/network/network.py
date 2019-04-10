import StarKillerMicrophysics as SKM
from StarKiller.interfaces import BurnType, EosType
from StarKiller.eos import Eos

class Network(object):
    def __init__(self):
        self.NetworkModule = SKM.Network()
        self.ActualNetworkModule = SKM.Actual_Network()
        self.RHSModule = SKM.actual_rhs_module
        self.eos = Eos()

        self.nspec = self.ActualNetworkModule.nspec
        self.nspec_evolve = self.ActualNetworkModule.nspec_evolve

        # These are python zero based indexes
        self.net_itemp = self.nspec_evolve
        self.net_ienuc = self.nspec_evolve + 1

        self.short_species_names = [self.NetworkModule.get_network_short_species_name(i+1).decode("ASCII").strip().lower() for i in range(self.nspec)]
        self.species_names = [self.NetworkModule.get_network_species_name(i+1).decode("ASCII").strip().lower() for i in range(self.nspec)]
        self.species_map = dict([(name, i) for i, name in enumerate(self.short_species_names)])

    def rhs(self, burn_state):
        eos_state = burn_state.to_eos_type()
        self.eos.evaluate(eos_state.eos_input_rt, eos_state)
        burn_state.from_eos_type(eos_state)
        self.RHSModule.actual_rhs(burn_state.state)

    def jacobian(self, burn_state):
        eos_state = burn_state.to_eos_type()
        self.eos.evaluate(eos_state.eos_input_rt, eos_state)
        burn_state.from_eos_type(eos_state)
        self.RHSModule.actual_jac(burn_state.state)
