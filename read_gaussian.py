import re
import numpy as np
import data_objects

def read_gaussian(file_name):
    lines = list(open(file_name))
    
    total_energy = 0.
    transition_energies = []
    transition_dipoles = []
    
    for enum, line in enumerate(lines):
        excitation_match = re.findall(r'Excited State   \d:', line)
        total_energy_match = re.findall(r'SCF Done:  E', line)
        
        if len(excitation_match) != 0:
            values = re.findall(r'\d+.\d+', line)
            transition_energies.append(float(values[0]) / 27.2114) #eV to hartree

        if len(total_energy_match) != 0:
            values = re.findall(r'-?\d+.\d+', line)
            total_energy = float(values[0])

        if "transition electric dipole moments (Au):" in line:
            for trns_dip in lines[enum+2:enum+7]:
                values = re.findall(r'-?\d+.\d+', trns_dip)
                transition_dipoles.append(np.array([float(x) for x in values[:3]]))
                
    return data_objects.DimerResult(total_energy, 
                                    transition_energies, 
                                    transition_dipoles)