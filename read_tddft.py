import re
import numpy as np
import data_objects

def read_monomer_TDDFT_file(file_name):    
    output_file = list(open(file_name))
                
    energy = None
    tdm = None
    is_here = True
    
    excited_state_lines = []

    for enum, line in enumerate(output_file):
        if "Excited State" in line or "Transition Dipole Moments (a.u.)" in line or "Oscillator Strengths (a.u.):" in line:
            excited_state_lines.append(enum)

    if(len(excited_state_lines) == 0):
        return (False, None, None)
            
    homo_lumo_state = 0
    homo_lumo_coeff = 0
    homo_lumo_state_line = None
    homo_lumo_dipole_line = None

    for i in range(len(excited_state_lines)-2):
        for line in output_file[excited_state_lines[i]:excited_state_lines[i+1]]:
            if "->" in line:
                MOs = re.findall(r'\d+', line)
                if abs(int(MOs[1]) - int(MOs[0])) != 1:
                    continue
                coeff = re.findall(r'-?\d+\.\d+', line)[0]
                if abs(float(coeff)) > homo_lumo_coeff:
                    homo_lumo_state = i+1
                    homo_lumo_coeff = float(coeff)
                    homo_lumo_state_line = output_file[excited_state_lines[i]]

    for line in output_file[excited_state_lines[-2]:excited_state_lines[-1]]:
        if re.match(fr' {homo_lumo_state}', line):
            homo_lumo_dipole_line = line

    nums = re.findall(r'\d+\.\d+', homo_lumo_state_line)
    energy= float(nums[0])
    tdm = np.array([float(x) for x in re.findall(r'-?\d+\.\d+', homo_lumo_dipole_line)])

        
    is_here = all([energy is not None, tdm is not None])
    
    if not is_here:
        print(file_name)
        return(is_here, None, None)
    else:
        return(is_here, energy, tdm)
    
    
def read_dimer_TDDFT_file(file_name):    
    output_file = list(open(file_name))
                
    total_energy = 0
    
    excited_state_lines = []
    transition_dipole_start = None
    oscillator_strength_start = None

    for enum, line in enumerate(output_file):
        if "Total ENERGY" in line:
            total_energy = re.findall(r'-?\d+\.\d+', line)
        
        if "Excited State" in line:
            excited_state_lines.append(enum)
            
        elif "Transition Dipole Moments (a.u.)" in line:
            transition_dipole_start = enum
            
        elif "Oscillator Strengths (a.u.):" in line:
            oscillator_strength_start = enum

    energies=[]
    tdms=[]
    
    for enum, line in enumerate(excited_state_lines):
        energies_strs = re.findall(r'-?\d*\.\d*', output_file[line])
        energies.append(float(energies_strs[0]))
        
        tdms_str = re.findall(r'-?\d*\.\d*', output_file[transition_dipole_start+enum+2])
        tdms.append(np.array([float(x) for x in tdms_str]))

    return data_objects.DimerResult(total_energy, energies, tdms)