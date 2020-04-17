import StarKillerMicrophysics as SKM
from StarKiller.interfaces import EosType

class BurnType(object):
    def __init__(self):
        self.BurnTypeModule = SKM.Burn_Type_Module()
        self.neqs = self.BurnTypeModule.neqs
        self.state = self.BurnTypeModule.burn_t()

    def to_eos_type(self):
        # Return an EosType object from this BurnType
        eos_state = EosType()
        self.BurnTypeModule.burn_to_eos(self.state, eos_state.state)
        return eos_state
        
    def from_eos_type(self, eos_state):
        # Given an EosType object, set this BurnType data
        self.BurnTypeModule.eos_to_burn(eos_state.state, self.state)

    def copy(self):
        # Return a deep copy of this object
        state_copy = BurnType()
        self.BurnTypeModule.copy_burn_t(state_copy.state, self.state)
        return state_copy

