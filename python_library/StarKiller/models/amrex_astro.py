import numpy as np
from StarKiller.network import Network

class AmrexAstroModel(object):
    def __init__(self, input_file):
        self.filename = input_file
        self.variables = []
        self.model_data = {}

        self.network = Network()
        print(self.network.species_names)

        self.filename = input_file
        self.read(self.filename)

    @property
    def fields(self):
        return list(self.model_data.keys())

    def data(self, field):
        try:
            assert(field in self.model_data.keys())
        except:
            try:
                field = self.network.shorten_species(field)
            except:
                raise
            else:
                pass
        return self.model_data[field]

    def reset(self):
        self.variables = ['radius']
        self.model_data = {}
        self.model_data['radius'] = []

    def read(self, input_file):
        self.reset()

        f = open(input_file, 'r')
        
        num_points_line = f.readline()
        num_points = int(num_points_line.split('=')[-1].strip())

        num_variables_line = f.readline()
        # add 1 for the radius variable
        num_variables = int(num_variables_line.split('=')[-1].strip()) + 1

        num_varnames_read = 0
        for l in f:
            ls = l.strip()

            if not ls:
                break
            elif ls[0] == '#':
                if num_varnames_read < 1:
                    num_varnames_read = 1
                else:
                    num_varnames_read += 1
                    
                # Read variable name if we're in the variable names section
                variable_name = ls[1:].strip()

                # Shorten species names to their abbreviations
                if variable_name in self.network.species_names:
                    variable_name = self.network.shorten_species(variable_name)
                
                self.variables.append(variable_name)
                self.model_data[variable_name] = []

                # Break if we have read all the variable names (radius not included)
                if num_varnames_read == num_variables - 1:
                    break

        fdata_entries = []
        for l in f:
            fdata_entries += l.strip().split()

        # Close file
        f.close()
        
        for ipt in range(num_points):
            for ivar in range(num_variables):
                ientry = ivar + ipt * (num_variables)
                variable_name = self.variables[ivar]
                self.model_data[variable_name].append(float(fdata_entries[ientry]))

        # Convert data to numpy arrays
        for vi in self.model_data.keys():
            self.model_data[vi] = np.array(self.model_data[vi])

