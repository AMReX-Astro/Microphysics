import StarKillerMicrophysics as SKM
from StarKiller.interfaces import EosType

BurnTypeModule = SKM.Burn_Type_Module()

class BurnType(object):
    neqs = BurnTypeModule.neqs
    
    def __init__(self):
        self.state = BurnTypeModule.burn_t()

    def to_eos_type(self):
        # Return an EosType object from this BurnType
        eos_state = EosType()
        BurnTypeModule.burn_to_eos(self.state, eos_state.state)
        return eos_state
        
    def from_eos_type(self, eos_state):
        # Given an EosType object, set this BurnType data
        BurnTypeModule.eos_to_burn(eos_state.state, self.state)

    def copy(self):
        # Return a deep copy of this object
        state_copy = BurnType()
        BurnTypeModule.copy_burn_t(state_copy.state, self.state)
        return state_copy

