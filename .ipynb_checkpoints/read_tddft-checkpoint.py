import re
import numpy as np
import data_objects

def read_dipole(lines):
    result = []
    
    for line in lines:
        result.append(float(re.findall(r'-?\d+\.\d+', line)[0]))
        
    return np.array(result)

def read_charges(lines):
    result = []
    
    for line in lines:
        result.append(float(re.findall(r'-?\d+\.\d+', line)[0]))
        
    return np.array(result)

def read_monomer_TDDFT_file(file_name):    
    output_file = list(open(file_name))
    
    n_atoms = 79
    
    total_energy = 0
    excitation_energy = 0
    
    ground_molecular_dipole = None
    transition_dipole = None
    excited_state_dipole = None
    
    ground_mulliken_charges = None
    ground_lowdin_charges = None
    
    transition_mulliken_charges = None
    transition_lowdin_charges = None
    
    excited_mulliken_charges = None
    excited_lowdin_charges = None
    
    for enum, line in enumerate(output_file):
        if "ground.energy" in line:
            total_energy = float(re.findall(r'-?\d+\.\d+', line)[0])

        elif "Excited State   1" in line:
            excitation_energy = float(re.findall(r'-?\d+\.\d+', line)[0])
            
        elif "ground.dipole" in line:
            ground_molecular_dipole = read_dipole(output_file[enum+2:enum+5])
            
        elif "excited_mulliken.transition_dipoles" in line:
            transition_dipole = read_dipole(output_file[enum+2:enum+5])
            
        elif "excited_mulliken.excited_state_dipoles" in line:
            excited_state_dipole = read_dipole(output_file[enum+2:enum+5])
            
        elif "excited_lowdin.transition_dipoles" in line:
            if((read_dipole(output_file[enum+2:enum+5]) != transition_dipole).all()):
                print(read_dipole(output_file[enum+2:enum+5]))
                print(transition_dipole)
            
        elif "excited_lowdin.excited_state_dipoles" in line:
            if((read_dipole(output_file[enum+2:enum+5]) != excited_state_dipole).all()):
                print(read_dipole(output_file[enum+2:enum+5]))
                print(excited_state_dipole)
            
        elif "ground_mulliken.charges" in line:
            ground_mulliken_charges = read_charges(output_file[enum+2:enum+2+n_atoms])
            
        elif "ground_lowdin.charges" in line:
            ground_lowdin_charges = read_charges(output_file[enum+2:enum+2+n_atoms])
            
        elif "excited_mulliken.transition_charges" in line:
            transition_mulliken_charges = read_charges(output_file[enum+2:enum+2+n_atoms])
            
        elif "excited_lowdin.transition_charges" in line:
            transition_lowdin_charges = read_charges(output_file[enum+2:enum+2+n_atoms])
                
        elif "excited_mulliken.excited_state_1_charges" in line:
            excited_mulliken_charges = read_charges(output_file[enum+2:enum+2+n_atoms])
                
        elif "excited_lowdin.excited_state_1_charges" in line:
            excited_lowdin_charges = read_charges(output_file[enum+2:enum+2+n_atoms])
    
    return data_objects.MonomerResult(total_energy = total_energy,
                                      molecular_dipole =  ground_molecular_dipole,
                                      mulliken_partial_charges = ground_mulliken_charges,
                                      lowdin_partial_charges = ground_lowdin_charges,
                                      
                                      transition_energy = excitation_energy,
                                      transition_dipole = transition_dipole,
                                      mulliken_transition_charges = transition_mulliken_charges,
                                      lowdin_transition_charges = transition_lowdin_charges,
                                      
                                      excited_molecular_dipole = excited_state_dipole,
                                      excited_mulliken_partial_charges = excited_mulliken_charges,
                                      excited_lowdin_partial_charges = excited_lowdin_charges)
    
def read_dimer_TDDFT_file(file_name):    
    output_file = list(open(file_name))
                
    total_energy = 0
    
    transition_energies = []
    transition_dipoles = []
    
    for enum, line in enumerate(output_file):
        if "single.energy" in line:
            total_energy = float(re.findall(r'-?\d+\.\d+', line)[0])
        
        elif "Excited State" in line:
            transition_energies.append(float(re.findall(r'-?\d+\.\d+', line)[0]))
            
        elif "Transition Dipole Moments (a.u.)" in line:
            for tdm_line in output_file[enum+2:enum+2+3]:
                tdms_str = re.findall(r'-?\d+\.\d+', tdm_line)
                transition_dipoles.append(np.array([float(x) for x in tdms_str]))

    return data_objects.DimerResult(total_energy, 
                                    transition_energies, 
                                    transition_dipoles)