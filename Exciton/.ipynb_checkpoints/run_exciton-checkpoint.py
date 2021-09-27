import subprocess
import os
import json
import numpy as np
import sys

from concurrent.futures import ProcessPoolExecutor
from simtk.openmm.app import *
from simtk.openmm import *

sys.path.insert(0, '..')
import read_pdbs


QCORE_PATH = os.environ["QCORE_PATH"]

def get_Mg_xyz(res, positions):
    for enum, atom in enumerate(res.atoms()):
        if atom.element.symbol == "Mg":
            xyz = positions[atom.index]
            
            return np.array([xyz[0]._value, xyz[1]._value, xyz[2]._value])

def Mg_Mg_distance(BCL_res1, BCL_res2, positions):
    positions = positions.in_units_of(unit.bohr)
    
    Mg1 = get_Mg_xyz(BCL_res1, positions)
    Mg2 = get_Mg_xyz(BCL_res2, positions)
    
    return np.linalg.norm(Mg1 - Mg2)

def read_hex_data(result, results_name, variable):
    
    hex_data = result[results_name][variable]["data"]["bytes"]
    byte_data = bytes.fromhex(" ".join([format(n, "02x") for n in hex_data]))
                                       
    dtype = result[results_name][variable]["dtype"]
    
    return np.frombuffer(byte_data, dtype=dtype)

def read_indexed_values(result, results_name, variable):
    
    data = result[results_name][variable]
    
    return data[0]
    

def run_monomer(qcore_string):
    """
    runs qcore with the input string
    """

    qcore_path = os.environ["QCORE_PATH"]

    options_str = " -n 1 -f json -s "
    
    job_str = qcore_path + options_str + "\"" + qcore_string + "\""

    json_run = subprocess.run(job_str,
                              shell=True,
                              stdout=subprocess.PIPE,
                              executable="/bin/bash",
                              universal_newlines=True)

    json_result = json.loads(json_run.stdout)
   
    total_energy = json_result["res"]["energy"]
    transition = json_result["res"]["excitation_1_energy"]
        
    return {
        "total_energy" : total_energy,
        "transition" : transition
    }
    
def run_dimer(qcore_string):
    """
    runs qcore with the input string
    """

    qcore_path = os.environ["QCORE_PATH"]

    options_str = " -n 1 -f json -s "
    
    job_str = qcore_path + options_str + "\"" + qcore_string + "\""

    json_run = subprocess.run(job_str,
                              shell=True,
                              stdout=subprocess.PIPE,
                              executable="/bin/bash",
                              universal_newlines=True)

    json_result = json.loads(json_run.stdout)
   
    eigenvalues = read_hex_data(json_result, "res", "eigenvalues")
    charge_centres = read_hex_data(json_result, "res", "charge_centres")

    couplings = read_indexed_values(json_result, "res", "couplings_inv_cm")
    distances = read_indexed_values(json_result, "res", "distances")
    
    eigenvectors = read_hex_data(json_result, "res", "eigenvectors")

    return {
        "eigenvalues" : eigenvalues.tolist(),
        "charge_centres" : charge_centres,
        "couplings" : couplings,
        "distances" : distances,
        "eigenvectors" : eigenvectors
    }
    
def set_dimer_data():
    return {
            "monomer A" : [],
            "monomer B" : [],
            "ring A" : [],
            "ring B" : [],
            "frame" : [],
            "distance" : [],
            "charge centre distance" : [],
            "coupling" : [],
            "exciton states" : [],
            "exciton transitions" : [],
            "eigenvectors" : []
        }
    
def set_monomer_data():
    return {
            "monomer" : [],
            "ring" : [],
            "frame" : [],
            "xtb states" : [],
            "xtb transition" : []
        }

def add_to_data(data, variable, value):
    data[variable].append(value)
    
    return data

def dump_data(data, dimer=True):
    if dimer:
        with open("exciton_dimer_data.json", 'w') as f:
            return json.dump(data, f)
    else:
        with open("exciton_monomer_data.json", 'w') as f:
            return json.dump(data, f)
    
def is_and_js():
    result = []
    
    for i in range(1, 28):
        for j in range(i+1, 28):
            result.append([i,j])
    
    return result
    
    
def write_qcore_str(xyzA, xyzB = None):
    if xyzB:
        qcore_str_template = "res := excitons(structure(xyz = {monomerA}) structure(xyz = {monomerB}) use_xtb = true hamiltonian = states)"
        return qcore_str_template.format(monomerA = xyzA, monomerB = xyzB)
    
    else:
        qcore_str_template = "res := xtb(structure(xyz = {monomerA}) model='chlorophyll')"
        return qcore_str_template.format(monomerA = xyzA)
        
    
if __name__ == "__main__":
    monomer_data = set_monomer_data()
    dimer_data = set_dimer_data()
    
    ring_assingments = json.load(open("../ring_assignment.json"))
    assign_ring = lambda i : ring_assingments["rings"][f"{i}"]
                
    monomer_indices = list(range(1, 28))
    dimer_indices = is_and_js()
    
    frame_range = range(1, 751, 250)
  

    for frame in frame_range:
        BCL_residues, positions = read_pdbs.read_pdb(f"../clean_pdbs/clean_md1_frame_{frame}.pdb")
            
        distances = list(map(lambda i_j : Mg_Mg_distance(BCL_residues[i_j[0]-1], BCL_residues[i_j[1]-1], positions), dimer_indices))
                    
        qcore_strings = map(lambda i_j : write_qcore_str(xyzA = read_pdbs.get_qcore_xyz(BCL_residues[i_j[0]-1], positions), xyzB = read_pdbs.get_qcore_xyz(BCL_residues[i_j[1]-1], positions)), dimer_indices)
            
        with ProcessPoolExecutor(max_workers=20) as pool:
            qcore_results = list(pool.map(run_dimer, qcore_strings))
    
        for enum, qcore_res in enumerate(qcore_results):
            i, j = dimer_indices[enum]

            row = [i, j, 
                   assign_ring(i), assign_ring(j),
                   frame,
                   distances[enum], qcore_res["distances"][1],
                   abs(qcore_res["couplings"][1]),
                   qcore_res["eigenvalues"],
                   [x - qcore_res["eigenvalues"][0] for x in qcore_res["eigenvalues"][1:]],
                   qcore-res["eigenvectors"]
                  ]
        
            for col, value in zip(dimer_data.keys(), row):
                dimer_data = add_to_data(dimer_data, col, value)
                
    dump_data(dimer_data)

    
    for frame in frame_range:
        BCL_residues, positions = read_pdbs.read_pdb(f"../clean_pdbs/clean_md1_frame_{frame}.pdb")
        
        qcore_strings = map(lambda i : write_qcore_str(xyzA = read_pdbs.get_qcore_xyz(BCL_residues[i-1], positions)), monomer_indices)
        
        with ProcessPoolExecutor(max_workers=20) as pool:
            qcore_results = list(pool.map(run_monomer, qcore_strings))
            
        for enum, qcore_res in enumerate(qcore_results):
            i = monomer_indices[enum]
            
            transition = qcore_res["transition"]
            
            row = [i, assign_ring(i), frame, [qcore_res["total_energy"], qcore_res["total_energy"] + transition], transition]
            
            for col, value in zip(monomer_data.keys(), row):
                monomer_data = add_to_data(monomer_data, col, value)
    
    dump_data(monomer_data, dimer=False)
    
    quit(0)
                
            
                
                
                
                
                
                
