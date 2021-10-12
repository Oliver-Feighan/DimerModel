import subprocess
import os
import json
import numpy as np
import sys

from concurrent.futures import ProcessPoolExecutor
from openmm.app import *
from openmm import *

sys.path.insert(0, '..')
sys.path.insert(0, '../Exciton')
import read_pdbs
import read_gaussian
import run_exciton

def write_qcore_str(xyzA, xyzB, partial_A, transition_A, excited_partial_A, partial_B, transition_B, excited_partial_B,
                    total_energy_A, total_energy_B, transition_energy_A, transition_energy_B):
    if xyzB:
        qcore_str_template = "res := excitons(structure(xyz = {monomerA}) structure(xyz = {monomerB}) coupling(charges = [['partial_charges', 1, {partial_A}], ['transition_charges', 1, {transition_A}], ['excited_partial_charges', 1, {excited_partial_A}], ['partial_charges', 2, {partial_B}], ['transition_charges', 2, {transition_B}], ['excited_partial_charges', 2, {excited_partial_B}]] method = charges) site_energies = [['total_energy', {total_energy_A}], ['total_energy', {total_energy_B}], ['transition_energy', {transition_energy_A}], ['transition_energy', {transition_energy_B}]] hamiltonian = states)"
        
        return qcore_str_template.format(monomerA = xyzA, monomerB = xyzB, partial_A = list(partial_A), transition_A = list(transition_A), excited_partial_A = list(excited_partial_A), partial_B = list(partial_B), transition_B = list(transition_B), excited_partial_B = list(excited_partial_B), total_energy_A = total_energy_A, total_energy_B = total_energy_B, transition_energy_A = transition_energy_A, transition_energy_B = transition_energy_B)
    
    else:
        qcore_str_template = "res := xtb(structure(xyz = {monomerA}) model='chlorophyll')"
        return qcore_str_template.format(monomerA = xyzA)

if __name__ == '__main__':
    ring_assingments = json.load(open("../ring_assignment.json"))
    assign_ring = lambda i : ring_assingments["rings"][f"{i}"]
    
    ij_indices = run_exciton.is_and_js()
    
    dimer_data = run_exciton.set_dimer_data()
    
    frame_range = range(1, 751, 250)
    
    for frame in frame_range:
        all_partial_charges = []
        all_transition_charges = []
        all_excited_partial_charges = []
    
        all_total_energies = []
        all_transition_energies = []
        
        for i in range(27):
            positions1, partial_charges = read_gaussian.read_positions_and_charges(f"../CAMB3LYP/CAMB3LYP_ground_trunc_bchla_{i+1}_frame_{frame}.log", 79, False)
            positions2, transition_charges = read_gaussian.read_positions_and_charges(f"../CAMB3LYP/CAMB3LYP_transition_trunc_bchla_{i+1}_frame_{frame}.log", 79, True)
            positions3, excited_partial_charges = read_gaussian.read_positions_and_charges(f"../CAMB3LYP/CAMB3LYP_excited_trunc_bchla_{i+1}_frame_{frame}.log", 79, False)
            
            tddft_result = read_gaussian.read_gaussian(f"../CAMB3LYP/CAMB3LYP_ground_trunc_bchla_{i+1}_frame_{frame}.log")
            
            all_partial_charges.append(partial_charges)
            all_transition_charges.append(transition_charges)
            all_excited_partial_charges.append(excited_partial_charges)
            all_total_energies.append(tddft_result.total_energy)
            all_transition_energies.append(tddft_result.excitations[0])

        BCL_residues, positions = read_pdbs.read_pdb(f"../clean_pdbs/clean_md1_frame_{frame}.pdb")
        
        distances = list(map(lambda i_j : run_exciton.Mg_Mg_distance(BCL_residues[i_j[0]-1], BCL_residues[i_j[1]-1], positions), ij_indices))

        
        qcore_strings = map(lambda i_j : write_qcore_str(xyzA = read_pdbs.get_qcore_xyz(BCL_residues[i_j[0]-1], positions), 
                                                         xyzB = read_pdbs.get_qcore_xyz(BCL_residues[i_j[1]-1], positions),
                                                         partial_A = all_partial_charges[i_j[0]-1],
                                                         transition_A = all_transition_charges[i_j[0]-1],
                                                         excited_partial_A = all_excited_partial_charges[i_j[0]-1],
                                                         partial_B = all_partial_charges[i_j[1]-1],
                                                         transition_B = all_transition_charges[i_j[1]-1],
                                                         excited_partial_B = all_excited_partial_charges[i_j[1]-1],
                                                         total_energy_A = all_total_energies[i_j[0]-1],
                                                         total_energy_B = all_total_energies[i_j[1]-1],
                                                         transition_energy_A = all_transition_energies[i_j[0]-1],
                                                         transition_energy_B = all_transition_energies[i_j[1]-1]), 
                            ij_indices)
        
        with ProcessPoolExecutor(max_workers=20) as pool:
            qcore_results = list(pool.map(run_exciton.run_dimer, qcore_strings))
    
        for enum, qcore_res in enumerate(qcore_results):
            i, j = ij_indices[enum]

            row = [i, j, 
                   assign_ring(i), assign_ring(j),
                   frame,
                   distances[enum], qcore_res["distances"][1],
                   abs(qcore_res["couplings"][1]),
                   qcore_res["eigenvalues"],
                   [x - qcore_res["eigenvalues"][0] for x in qcore_res["eigenvalues"][1:]],
                   qcore_res["eigenvectors"]
                  ]
        
            for col, value in zip(dimer_data.keys(), row):
                dimer_data = run_exciton.add_to_data(dimer_data, col, value)
                
    run_exciton.dump_data(dimer_data, "exciton_dimer_data.json")
        