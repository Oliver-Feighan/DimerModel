import re
import numpy as np
import data_objects

def read_positions_and_charges(file, n_atoms, is_transition):
    nuclear_charge = {
        "Mg" : 12,
        "C" : 6,
        "N" : 7,
        "O" : 8,
        "H" : 1
    }

    positions = []
    charges = []

    with open(file) as f:
        lines = list(f.readlines())
        for enum, line in enumerate(lines):
            if " Symbolic Z-matrix:" in line:
                for position_line in lines[enum+2:enum+2+n_atoms]:
                    element = re.findall(r'[a-zA-Z]+', position_line)[0]
                    coords = np.array([float(x) for x in re.findall(r'-?\d+\.\d+|-?\d+\.', position_line)])
                    if len(coords) != 3:
                        print(coords, position_line)

                    positions.append(coords)

            if "Mulliken charges:" in line:
                for charge_line in lines[enum+2:enum+2+n_atoms]:
                    population = float(re.findall(r'-?\d+\.\d+', charge_line)[0])
                    element = re.findall(r'[a-zA-Z]+', charge_line)[0]
                    
                    if is_transition:
                        charges.append(population - nuclear_charge[element])
                    else:
                        charges.append(population)


    positions = np.array(positions)
    charges = np.array(charges)
    
    return positions, charges

def get_excited_state_lines(lines):
    result = []
    
    for enum, line in enumerate(lines):
        if re.search(r'Excited State   \d:', line) or re.search(r'SavETr', line):
            result.append(enum)
            
    return result

def read_gaussian(file_name):
    lines = list(open(file_name))
    
    total_energy = 0.
    transition_energies = []
    transition_characters = []
    transition_dipoles = []
    
    excited_state_lines = get_excited_state_lines(lines)
    
    for i in range(len(excited_state_lines)-1):
        character = []
        for line in lines[excited_state_lines[i]:excited_state_lines[i+1]]:
            if re.search(r'->', line) or re.search(r'<-', line):
                character.append([float(x) for x in re.findall(r'-?\d+\.?\d+', line)])
                
        transition_characters.append(character)
                    
    
    for enum, line in enumerate(lines):
        total_energy_match = re.findall(r'SCF Done:  E', line)
        
        if enum in excited_state_lines:
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
                                    transition_dipoles,
                                    transition_characters)