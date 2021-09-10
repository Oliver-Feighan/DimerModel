from simtk.openmm.app import *
from simtk.openmm import *

import glob
import re

def get_BCL_residues(topology):
    result = []
    
    for residue in topology.residues():
        if residue.name == "BCL":
            result.append(residue)

    return result    
    
def extract_xyz_line(atom, positions, convert_symbol=False):
    symbol = atom.element.symbol
    x = positions[atom.index][0]._value
    y = positions[atom.index][1]._value
    z = positions[atom.index][2]._value
    
    if not convert_symbol:
        return f"{symbol:2s} {x:3.5f} {y:3.5f} {z:3.5f}\n"
    else: 
        return f"{convert_symbol:2s} {x:3.5f} {y:3.5f} {z:3.5f}\n"
    
def write_monomer(residue, positions, file_name, truncate=True):
    positions = positions.in_units_of(unit.angstrom)
    
    index_range = lambda x, y : list(range(x, y + 1))
    
    tail_indices = index_range(13, 15) + index_range(46, 65) + index_range(101, 139)
    C_to_H_index = 13
    
    xyz = ""
    n_atoms = 0
    
    for enum, atom in enumerate(residue.atoms()):
        if enum not in tail_indices:
            xyz += extract_xyz_line(atom, positions)
            n_atoms += 1
            
        if enum == C_to_H_index:
            xyz += extract_xyz_line(atom, positions, "H")
            n_atoms += 1
            
    header = f"{n_atoms}\n\n"
    
    file = open(f"monomer_xyzs/{file_name}", 'w')
    print(header+xyz, file=file)
    file.close()

def write_dimer(residue1, residue2, positions, file_name, truncate=True):
    positions = positions.in_units_of(unit.angstrom)
    
    index_range = lambda x, y : list(range(x, y + 1))
    
    tail_indices = index_range(13, 15) + index_range(46, 65) + index_range(101, 139)
    C_to_H_index = 13
    
    xyz = ""
    n_atoms = 0
    
    for enum, atom in enumerate(residue1.atoms()):
        if enum not in tail_indices:
            xyz += extract_xyz_line(atom, positions)
            n_atoms += 1
            
        if enum == C_to_H_index:
            xyz += extract_xyz_line(atom, positions, "H")
            n_atoms += 1
    
    for enum, atom in enumerate(residue2.atoms()):
        if enum not in tail_indices:
            xyz += extract_xyz_line(atom, positions)
            n_atoms += 1
            
        if enum == C_to_H_index:
            xyz += extract_xyz_line(atom, positions, "H")
            n_atoms += 1
    
    header = f"{n_atoms}\n\n"
    
    file = open(f"dimer_xyzs/{file_name}", 'w')
    print(header+xyz, file=file)
    file.close()

    
if __name__ == "__main__":
    n_frames = 500
    
    for i in range(n_frames):
        frame_number = list(range(1, 50001, 50))[int(1000/n_frames) * i]
    
        print(frame_number)
    
        file = f"clean_pdbs/clean_md1_frame_{frame_number}.pdb"
        
        pdb = PDBFile(file)

        positions = pdb.getPositions(asNumpy = True)
        topology = pdb.getTopology()

        BCL_residues = get_BCL_residues(topology)

        for enum1, residue1 in enumerate(BCL_residues):
            write_monomer(residue1, positions, f"trunc_bchla_{enum1+1}_frame_{frame_number}.xyz")
            
            for enum2, residue2 in enumerate(BCL_residues):
                if enum1 >= enum2:
                    continue
                
                write_dimer(residue1, residue2, positions, f"trunc_bchla_{enum1+1}_bchla_{enum2+1}_frame_{frame_number}.xyz")