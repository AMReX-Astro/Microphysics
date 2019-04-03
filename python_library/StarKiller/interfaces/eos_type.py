import StarKillerMicrophysics as SKM

EosTypeModule = SKM.Eos_Type_Module()

class EosType(object):
    eos_input_rt = EosTypeModule.eos_input_rt
    eos_input_rh = EosTypeModule.eos_input_rh
    eos_input_tp = EosTypeModule.eos_input_tp
    eos_input_rp = EosTypeModule.eos_input_rp
    eos_input_re = EosTypeModule.eos_input_re
    eos_input_ps = EosTypeModule.eos_input_ps
    eos_input_ph = EosTypeModule.eos_input_ph
    eos_input_th = EosTypeModule.eos_input_th    
    
    def __init__(self):
        self.state = EosTypeModule.eos_t()
