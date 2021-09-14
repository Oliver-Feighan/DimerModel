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
    

def run_qcore(qcore_string):
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
    
    return {
        "eigenvalues" : eigenvalues.tolist(),
        "charge_centres" : charge_centres,
        "couplings" : couplings,
        "distances" : distances
    }
    
    
def set_data():
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
    }

def add_to_data(data, variable, value):
    data[variable].append(value)
    
    return data

def dump_data(data):
    with open("exciton_data.json", 'w') as f:
        return json.dump(data, f)

    
def is_and_js():
    result = []
    
    for i in range(1, 28):
        for j in range(i+1, 28):
            result.append([i,j])
    
    return result
    
    
def write_qcore_str(xyzA, xyzB):
    qcore_str_template = "res := excitons(structure(xyz = {monomerA}) structure(xyz = {monomerB}) use_xtb = true hamiltonian = states)"

    return qcore_str_template.format(monomerA = xyzA, monomerB = xyzB)
    

    
if __name__ == "__main__":
    data = set_data()
    
    ring_assingments = json.load(open("../ring_assignment.json"))
    assign_ring = lambda i : ring_assingments["rings"][f"{i}"]
                
    index_tuples = is_and_js()
    
    for frame in range(1, 49951, 250):
        BCL_residues, positions = read_pdbs.read_pdb(f"../clean_pdbs/clean_md1_frame_{frame}.pdb")
            
        distances = list(map(lambda i_j : Mg_Mg_distance(BCL_residues[i_j[0]-1], BCL_residues[i_j[1]-1], positions), index_tuples))
                    
        qcore_strings = map(lambda i_j : write_qcore_str(xyzA = read_pdbs.get_qcore_xyz(BCL_residues[i_j[0]-1], positions), xyzB = read_pdbs.get_qcore_xyz(BCL_residues[i_j[1]-1], positions)), index_tuples)
            
        with ProcessPoolExecutor(max_workers=20) as pool:
            qcore_results = list(pool.map(run_qcore, qcore_strings))
    
        for enum, qcore_res in enumerate(qcore_results):
            i, j = index_tuples[enum]

            row = [i, j, assign_ring(i), assign_ring(j), frame, distances[enum], qcore_res["distances"][1], abs(qcore_res["couplings"][1]), qcore_res["eigenvalues"], [x - qcore_res["eigenvalues"][0] for x in qcore_res["eigenvalues"][1:]]]
        
            for col, value in zip(data.keys(), row):
                data = add_to_data(data, col, value)
                
    dump_data(data)
    
    quit(0)
                
            
                
                
                
                
                
                
