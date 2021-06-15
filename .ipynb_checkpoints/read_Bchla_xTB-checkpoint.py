import numpy as np
import re

def read_monomer_Bchla_xTB_file(file_name):
    output_file = list(open(file_name))
    
    energy = None
    tdm = None
    transition_charges = []
    is_here = None

    transition_charge_line = None
    
    homo_lumo = False
    
    for enum, line in enumerate(output_file):
        if "Excitation 0 1" in line:
            homo_lumo = True
        
        if "Excitation 0 2" in line:
            homo_lumo = False
    
        if not homo_lumo:
            continue
            
        if "Excitation energy:" in line:
            energy = float(re.findall(r'-?\d+\.\d+', line)[0])
        elif "Transition dipole" in line:
            tdm = np.array([float(x) for x in re.findall(r'-?\d+\.\d+', line)])
        elif "transition charges and positions" in line:
            transition_charge_line = enum
        
        if transition_charge_line is not None:
            if enum > transition_charge_line+3 and len(transition_charges) < 79:
                transition_charges.append(float(re.findall(r'-?\d+\.\d+', line)[0]))
                
         
            
    is_here = energy is not None and tdm is not None
            
    return (is_here, energy, tdm, transition_charges)

def read_dimer_Bchla_xTB_file(file_name):
    output_file = list(open(file_name))
    
    energies = []
    tdms = []
    is_here = None
    
    for enum, line in enumerate(output_file):
        if "Excitation energy:" in line:
            energies.append(float(re.findall(r'-?\d+\.\d+', line)[0]))
        elif "Transition dipole" in line:
            tdms.append(np.array([float(x) for x in re.findall(r'-?\d+\.\d+', line)]))
            
    is_here = energies and tdms
            
    return (is_here, energies, tdms)