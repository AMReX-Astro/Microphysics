import StarKillerMicrophysics as SKM
from StarKiller.interfaces import EosType

EosModule = SKM.Eos_Module()

class Eos(object):
    @staticmethod
    def evaluate(input_mode, input_state, use_raw_inputs=False):
        EosModule.eos(input_mode, input_state.state, use_raw_inputs)
