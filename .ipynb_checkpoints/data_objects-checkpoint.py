import re
import numpy as np

class Structure:
    Na_index = 14
    Nc_index = 31

    def __init__(self, bchla, frame, given_file_name = ""):
        self.symbols = []
        self.coords = []
        
        angstrom_to_bohr = 1.88973
                
        file_name = None
            
        if given_file_name != "":
            file_name = given_file_name
        else:
            file_name = f"monomer_xyzs/trunc_bchla_{bchla}_frame_{frame}.xyz"
        
        lines = list(open(file_name))
        lines = lines[2:]

        for line in lines:
            symbol = re.findall(r'[a-zA-Z]+', line)
            coord  = np.array([float(y) * angstrom_to_bohr for y in re.findall(r'-?\d+.\d+', line)])
            
            if len(symbol) == 0 or len(coord) == 0:
                continue
                
            self.symbols.append(symbol[0])
            self.coords.append(coord)
    
        assert(len(self.symbols) == 79)
        assert(len(self.coords) == 79)
    
        self.center = self.coords[0]

        self.Na_Nc = self.coords[self.Na_index] - self.coords[self.Nc_index]

class DimerResult():
    def __init__(self, total_energy, excitations, transition_dipoles):
        self.total_energy = total_energy
        
        self.excitations = excitations
        self.state_energies = [total_energy + x for x in self.excitations]
        
        self.transition_dipoles = transition_dipoles
        
class MonomerResult():
    
    def __init__(self, total_energy, molecular_dipole, partial_charges, transition_energy, transition_dipole, transition_charges, excited_molecular_dipole, excited_partial_charges):
        
        self.total_energy = total_energy
        self.molecular_dipole = molecular_dipole
        self.partial_charges = partial_charges

        self.transition_energy = transition_energy
        self.transition_dipole = transition_dipole
        self.transition_charges = transition_charges
        
        self.excited_molecular_dipole = excited_molecular_dipole
        self.excited_partial_charges = excited_partial_charges

        
        
class FrenkelResult():
    def __init__(self, hamiltonian, tdm_angle, eigvals, eigvecs):
        self.hamiltonian = hamiltonian
        self.tdm_angle = tdm_angle
        self.eigvals = eigvals
        self.eigvecs = eigvecs
        