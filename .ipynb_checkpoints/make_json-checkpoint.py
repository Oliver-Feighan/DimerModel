import json
import pickle as pkl
import glob
import numpy as np
import re
import matplotlib.pyplot as plt 

import read_tddft
import read_Bchla_xTB

class Structure:
    Na_index = 14
    Nc_index = 31

    def __init__(self, bchla, frame):
        self.symbols = []
        self.coords = []
        
        angstrom_to_bohr = 1.88973
        
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

def angle(vec1, vec2):
    num = np.dot(vec1, vec2)
    dom = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    
    angle = np.rad2deg(np.arccos(num/dom))
    
    return min(angle, 180-angle)

def dipole_coupling(tdm1, tdm2, bchla1_structure, bchla2_structure):
    r_vec = bchla1_structure.center - bchla2_structure.center
    r_mag = np.linalg.norm(r_vec)
    
    Nac1 = bchla1_structure.Na_Nc * (np.linalg.norm(tdm1)/np.linalg.norm(bchla1_structure.Na_Nc))
    Nac2 = bchla1_structure.Na_Nc * (np.linalg.norm(tdm2)/np.linalg.norm(bchla2_structure.Na_Nc))
    
    cos1 = np.dot(tdm1, r_vec) / (np.linalg.norm(tdm1) * r_mag)
    cos2 = np.dot(tdm2, r_vec) / (np.linalg.norm(tdm2) * r_mag)
    cos12 = np.dot(tdm1, tdm2) / (np.linalg.norm(tdm1) * np.linalg.norm(tdm2))
    
    vec_v = (np.dot(tdm1, tdm2)/(r_mag ** 3)) - 3*((np.dot(tdm1, r_vec) * np.dot(tdm2, r_vec))/(r_mag **5))
    
    return vec_v

def charge_coupling(tcs1, tcs2, bchla1_structure, bchla2_structure):
    result = 0
    
    for n_i, i in enumerate(bchla1_structure.coords):
        for n_j, j in enumerate(bchla2_structure.coords):
            r = np.linalg.norm(i-j)
            
            result += tcs1[n_i] * tcs2[n_j] / r
            
    return float(result)

def switch_coupling(tdm1, tdm2, tcs1, tcs2, bchla1_structure, bchla2_structure, method):
    if method == "dipole":
        return dipole_coupling(tdm1, tdm2, bchla1_structure, bchla2_structure)
    elif method == "charges":
        return charge_coupling(tcs1, tcs2, bchla1_structure, bchla2_structure)

def run_Frenkel_Hamiltonian(energy1, energy2, tdm1, tdm2, tcs1, tcs2, bchla1_structure, bchla2_structure, method):
    v = switch_coupling(tdm1, tdm2, tcs1, tcs2, bchla1_structure, bchla2_structure, method)
    tdms_angle = angle(tdm1, tdm2)

    mat = np.array([energy1, v, v, energy2])
    mat = mat.reshape((2,2))
        
    eigvals, eigvecs = np.linalg.eig(mat)

    return eigvals, v, tdms_angle
    
def make_PBE0_json():
    result = {}
    
    for file in glob.glob("PBE0/PBE0_*out"):
        nums = re.findall(r'\d+', file)
        frame = nums[-1]

        if f"frame_{frame}" not in result:
            result[f"frame_{frame}"] = {}
        
        n_mers = len(re.findall(r'bchla', file))
        
        if n_mers == 2:
            bchla1 = nums[-3]
            bchla2 = nums[-2]
            is_here_dimer, dimer_energies, dimer_tdms = read_tddft.read_dimer_TDDFT_file(file)
            is_here1, energy1, tdm1 = read_tddft.read_monomer_TDDFT_file(f"PBE0/PBE0_trunc_bchla_{bchla1}_frame_{frame}.out")
            is_here2, energy2, tdm2 = read_tddft.read_monomer_TDDFT_file(f"PBE0/PBE0_trunc_bchla_{bchla2}_frame_{frame}.out")
            
            if not is_here_dimer or not is_here1 or not is_here2:
                continue
            
            bchla1_structure = Structure(bchla1, frame)
            bchla2_structure = Structure(bchla2, frame)
            
            dipole_exciton_energies, dipole_coupling, dipole_angle = run_Frenkel_Hamiltonian(energy1, energy2, tdm1, tdm2, [], [], bchla1_structure, bchla2_structure, "dipole")

            #charges_exciton_energies, charges_coupling, dipole_angle = run_Frenkel_Hamiltonian(energy1, energy2, tdm1, tdm2, tcs1, tcs2, bchla1_structure, bchla2_structure, "charges")

            distance = np.linalg.norm(bchla1_structure.center - bchla2_structure.center)
            
            result[f"frame_{frame}"][f"bchla_{bchla1}_bchla_{bchla2}"] = {"dimer_energies" : dimer_energies,
                                                                          "dimer_tdms" : [x.tolist() for x in dimer_tdms],
                                                                          "dipole_exciton_energies" : dipole_exciton_energies.tolist(),
                                                                          "distance" : distance,
                                                                          "dipole_coupling" : dipole_coupling,
                                                                          "angle"    : dipole_angle,
                                                                          "monomer_energies" : [energy1, energy2]}
                

    return {"PBE0" : result}
    
def make_Bchla_xTB_json():
    result = {}
    
    for file in glob.glob("Bchla_xTB/Bchla_*out"):
        nums = re.findall(r'\d+', file)
        frame = nums[-1]

        if f"frame_{frame}" not in result:
            result[f"frame_{frame}"] = {}
        
        n_mers = len(re.findall(r'bchla', file))
        
        if n_mers == 2:
            bchla1 = nums[-3]
            bchla2 = nums[-2]
            is_here_dimer, dimer_energies, dimer_tdms = read_Bchla_xTB.read_dimer_Bchla_xTB_file(file)
            is_here1, energy1, tdm1, tcs1 = read_Bchla_xTB.read_monomer_Bchla_xTB_file(f"Bchla_xTB/Bchla_trunc_bchla_{bchla1}_frame_{frame}.out")
            is_here2, energy2, tdm2, tcs2 = read_Bchla_xTB.read_monomer_Bchla_xTB_file(f"Bchla_xTB/Bchla_trunc_bchla_{bchla2}_frame_{frame}.out")
            
            bchla1_structure = Structure(bchla1, frame)
            bchla2_structure = Structure(bchla2, frame)
            
            dipole_exciton_energies, dipole_coupling, dipole_angle = run_Frenkel_Hamiltonian(energy1, energy2, tdm1, tdm2, tcs1, tcs2, bchla1_structure, bchla2_structure, "dipole")

            charges_exciton_energies, charges_coupling, dipole_angle = run_Frenkel_Hamiltonian(energy1, energy2, tdm1, tdm2, tcs1, tcs2, bchla1_structure, bchla2_structure, "charges")

            
            distance = np.linalg.norm(bchla1_structure.center - bchla2_structure.center)
            
            result[f"frame_{frame}"][f"bchla_{bchla1}_bchla_{bchla2}"] = {"dimer_energies" : dimer_energies,
                                                                          "dimer_tdms" : [x.tolist() for x in dimer_tdms],
                                                                          "dipole_exciton_energies" : dipole_exciton_energies.tolist(),
                                                                          "charges_exciton_energies" : charges_exciton_energies.tolist(),
                                                                          "distance" : distance,
                                                                          "dipole_coupling" : dipole_coupling,
                                                                          "charges_coupling" : charges_coupling,
                                                                          "angle"    : dipole_angle,
                                                                          "monomer_energies" : [energy1, energy2]}
                

    return {"Bchla_xTB" : result}

def make_json(method):
    if method == "Bchla_xTB":
        return make_Bchla_xTB_json()
    elif method == "PBE0":
        return make_PBE0_json()


def extract_bchla_nums(bchlas):
    nums = re.findall(r'\d+', bchlas)
    return [int(x) for x in nums]

if __name__ == "__main__":
    #bchla_results = make_json("Bchla_xTB")
    
    #bchla_result_file = open("Bchla_xTB_results.json", 'w')
    #json.dump(bchla_results, bchla_result_file, indent=4)
    #bchla_result_file.close()
    
    PBE0_results = make_json("PBE0")
    
    PBE0_result_file = open("PBE0_results.json", 'w')
    json.dump(PBE0_results, PBE0_result_file, indent=4)
    PBE0_result_file.close()    
    
    quit()
    
    fig, ax = plt.subplots()
    
    counter = 0
    
    for frame in list(dimers.keys()):
        if frame in monomers:
            for bchlas in list(dimers[frame].keys()):
                bchla1, bchla2 = extract_bchla_nums(bchlas)
                
                dimer_energies = dimers[frame][bchlas]["energies"]
                bchla1_energy = monomers[frame][f"bchla_{bchla1}"]["energy"]
                bchla2_energy = monomers[frame][f"bchla_{bchla2}"]["energy"]
                
                if any([dimer_energies is None, bchla1_energy is None, bchla2_energy is None]):
                    continue
                    
                for i in dimer_energies:
                    ax.scatter(counter, i, color='black', marker='|')
                
                for i in [bchla1_energy, bchla2_energy]:
                    ax.scatter(counter, i, color='red', marker='x')
                
                counter += 1

    plt.show()