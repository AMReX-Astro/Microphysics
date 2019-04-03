import StarKillerMicrophysics as SKM
from StarKiller.interfaces import BurnType, EosType
from StarKiller.eos import Eos

NetworkModule = SKM.Network()
ActualNetworkModule = SKM.Actual_Network()
RHSModule = SKM.actual_rhs_module

class Network(object):
    nspec = ActualNetworkModule.nspec
    nspec_evolve = ActualNetworkModule.nspec_evolve

    # These are python zero based indexes
    net_itemp = nspec_evolve
    net_ienuc = nspec_evolve + 1
    
    short_species_names = [NetworkModule.get_network_short_species_name(i+1).decode("ASCII").strip().lower() for i in range(nspec)]
    species_names = [NetworkModule.get_network_species_name(i+1).decode("ASCII").strip().lower() for i in range(nspec)]

    @staticmethod
    def rhs(burn_state):
        Eos.evaluate(EosType.eos_input_rt, burn_state)
        RHSModule.actual_rhs(burn_state.state)

    @staticmethod
    def jacobian(burn_state):
        Eos.evaluate(EosType.eos_input_rt, burn_state)
        RHSModule.actual_jac(burn_state.state)
