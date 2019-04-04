import StarKillerMicrophysics as SKM
from StarKiller.interfaces import EosType

class Eos(object):
    def __init__(self):
        self.EosModule = SKM.Eos_Module()

    def evaluate(self, input_mode, input_state, use_raw_inputs=False):
        self.EosModule.eos(input_mode, input_state.state, use_raw_inputs)
