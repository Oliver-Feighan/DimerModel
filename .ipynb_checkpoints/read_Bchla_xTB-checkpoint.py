import numpy as np
import re

import data_objects

def read_monomer_Bchla_xTB_file(file_name):
    output_file = list(open(file_name))
        
    total_energy = None
    molecular_dipole = []
    partial_charges = []
    
    transition_energy = None
    transition_dipole = []
    transition_charges = []
    
    excited_molecular_dipole = []
    excited_partial_charges = []
    
    
    for enum, line in enumerate(output_file):
        if "res.energy" in line:
            total_energy = float(re.findall(r'-?\d+\.\d+', line)[0])
        
        elif "res.atomic_charges" in line:            
            for charge_line in output_file[enum+2:enum+2+79]: #2 for the lines until the charge list, 79 for the number of atoms
                partial_charges.append(float(re.findall(r'-?\d+\.\d+', charge_line)[0]))
        
        elif "res.dipole" in line:
            for dipole_line in output_file[enum+2:enum+2+3]:
                molecular_dipole.append(float(re.findall(r'-?\d+\.\d+', dipole_line)[0]))
        
        elif "res.excitation_1_energy" in line:
            transition_energy = float(re.findall(r'-?\d+\.\d+', line)[0])
            
        elif "res.excitation_1_transition_dipole" in line:
            for transition_dipole_line in output_file[enum+1:enum+1+3]:
                transition_dipole.append(float(re.findall(r'-?\d+\.\d+', transition_dipole_line)[0]))
        
        elif "res.excitation_1_transition_charges" in line:
            for transition_charges_line in output_file[enum+2:enum+2+79]:
                transition_charges.append(float(re.findall(r'-?\d+\.\d+', transition_charges_line)[0]))
                
        elif "res.excitation_1_excited_dipole" in line:
            for excited_dipole_line in output_file[enum+1:enum+1+3]:
                excited_molecular_dipole.append(float(re.findall(r'-?\d+\.\d+', excited_dipole_line)[0]))
        
        elif "res.excitation_1_excited_atomic_charges" in line:
            for excited_charges_line in output_file[enum+2:enum+2+79]:
                excited_partial_charges.append(float(re.findall(r'-?\d+\.\d+', excited_charges_line)[0]))
                    
    return data_objects.MonomerResult(total_energy, molecular_dipole, partial_charges, transition_energy, transition_dipole, transition_charges, excited_molecular_dipole, excited_partial_charges)
    
def read_dimer_Bchla_xTB_file(file_name):
    output_file = list(open(file_name))
    
    total_energy = None
    
    energies = []
    tdms = []
    is_here = None
    
    for enum, line in enumerate(output_file):
        if "res.energy" in line:
            total_energy = float(re.findall(r'-?\d+\.\d+', line)[0])
        
        if re.search("res.excitation_\d_energy", line):
            energies.append(float(re.findall(r'-?\d+\.\d+', line)[0]))
        elif re.search("res.excitation_\d_transition_dipole", line):
            tdm = []
            for transition_dipole_line in output_file[enum+1:enum+1+3]:
                tdm.append(float(re.findall(r'-?\d+\.\d+', transition_dipole_line)[0]))
                
            tdms.append(tdm)
            
    is_here = energies and tdms
            
    return data_objects.DimerResult(total_energy, energies, tdms)