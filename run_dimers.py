import numpy as np
import simtk.openmm.app
import simtk.openmm
import os
import subprocess
import json
import matplotlib.pyplot as plt

class single_result():
	def __init__(self, name, energy, transition_dipole, transition_charges):
		self.name = name
		self.energy = energy
		self.transition_dipole = transition_dipole
		self.transition_charges = transition_charges

def run_qcore(xyz_str):
	qcore_path = os.environ["qcore_path"]

	parameters = {'k_s': 1.602, 'k_p': 3.328, 'Mg_s': 1.152, 'Mg_p': 1.217, 'N_s': 1.088, 'N_p': 0.895, 'a_x': 0.053, 'y_K': 1.945, 'y_J': 3.955}

	qcore_str = f"{qcore_path} -n 1 -f json --schema none -s \" res := bchla(structure(xyz = {xyz_str}) model='gfn1' input_params={parameters})\""

	json_run = subprocess.run(qcore_str,
		shell=True,
		stdout=subprocess.PIPE,
		executable="/bin/bash",
		universal_newlines=True)

	result = json.loads(json_run.stdout)

	return single_result("res", result["res"]["excitation_energy"], result["res"]["transition_dipole"], result["res"]["transition_charges"])


def charges_coupling(residue1, residue2, positions, result1, result2):
	positions = positions.in_units_of(simtk.openmm.unit.bohr)
	result = 0

	for enum1, atom1 in enumerate(residue1.atoms()):
		for enum2, atom2 in enumerate(residue2.atoms()):
			pos1 = positions[atom1.index]
			pos2 = positions[atom2.index]

			r = np.linalg.norm(pos1._value - pos2._value)

			q1_q2 = result1.transition_charges[enum1] * result2.transition_charges[enum2]

			result += q1_q2/r

	return result

def get_distance(residue1, residue2, positions):
	positions = positions.in_units_of(simtk.openmm.unit.bohr)

	Mg1 = None

	Mg2 = None

	for atom in residue1.atoms():
		if atom.name == "MG":
			Mg1 = atom

	for atom in residue2.atoms():
		if atom.name == "MG":
			Mg2 = atom

	r_vec = positions[Mg1.index]._value - positions[Mg2.index]._value
	r_mag = np.linalg.norm(r_vec)

	return r_mag

def dipole_coupling(residue1, residue2, positions, result1, result2, orient=False):
	positions = positions.in_units_of(simtk.openmm.unit.bohr)
	
	tdm1 = result1.transition_dipole
	tdm2 = result2.transition_dipole

	Mg1 = None
	Na1 = None
	Nc1 = None

	Mg2 = None
	Na2 = None
	Nc2 = None

	for atom in residue1.atoms():
		if atom.name == "MG":
			Mg1 = atom
		elif atom.name == "N_A":
			Na1 = atom
		elif atom.name == "N_C":
			Nc1 = atom

	for atom in residue2.atoms():
		if atom.name == "MG":
			Mg2 = atom
		elif atom.name == "N_A":
			Na2 = atom
		elif atom.name == "N_C":
			Nc2 = atom

	r_vec = positions[Mg1.index]._value - positions[Mg2.index]._value
	r_mag = np.linalg.norm(r_vec)

	if orient:
		NaNc1 = positions[Na1.index]._value - positions[Nc1.index]._value
		NaNc1_unit = NaNc1 / np.linalg.norm(NaNc1)
		tdm1 = np.linalg.norm(tdm1) * NaNc1_unit

		NaNc2 = positions[Na2.index]._value - positions[Nc2.index]._value
		NaNc2_unit = NaNc2 / np.linalg.norm(NaNc2)
		tdm2 = np.linalg.norm(tdm2) * NaNc2_unit


	v = np.dot(tdm1, tdm2)/(r_mag ** 3) - 3 * ((np.dot(tdm1, r_vec) * np.dot(tdm2, r_vec))/ (r_mag ** 5))

	return v

def resolve_coupling(residue1, residue2, positions, result1, result2, coupling):
	if coupling == "charges":
		return charges_coupling(residue1, residue2, positions, result1, result2)
	
	elif coupling == "raw_dipole":
		return dipole_coupling(residue1, residue2, positions, result1, result2)

	elif coupling == "NaNc_dipole":
		return dipole_coupling(residue1, residue2, positions, result1, result2, True)


def frenkel_hamiltonian_energy(residue1, residue2, positions, result1, result2, coupling):
	e1 = result1.energy
	e2 = result2.energy
	v = resolve_coupling(residue1, residue2, positions, result1, result2, coupling)

	mat = np.array([e1, v, v, e2])
	mat = mat.reshape((2,2))

	eigvals, eigvecs = np.linalg.eig(mat)
	return min(eigvals)

def read_monomer_result(enum):
	qcore_result = json.load(open(f"monomers/Bchla_{enum+1}.json"))

	return single_result(f"Bchla_{enum+1}", qcore_result["res"]["excitation_energy"], qcore_result["res"]["transition_dipole"], qcore_result["res"]["transition_charges"])

def read_dimer_result(enum1, enum2):
	qcore_result = json.load(open(f"dimers/Bchla_{enum1+1}_Bchla_{enum2+1}.json"))

	return single_result(f"Bchla_{enum1+1}_Bchla_{enum2+1}", qcore_result["res"]["excitation_energy"], qcore_result["res"]["transition_dipole"], qcore_result["res"]["transition_charges"])


if __name__ == '__main__':
	pdb = simtk.openmm.app.pdbfile.PDBFile('clean_md1.pdb')

	topology = pdb.getTopology()
	positions = pdb.getPositions(asNumpy = True)
	positions = positions.in_units_of(simtk.openmm.unit.angstrom)

	bchla_residues = []
	bchla_xyzs = []
	bchla_results = []

	distances = []
	dimer_energies = []
	charge_coupling_energies = []
	dipole_coupling_energies = []
	orient_coupling_energies = []



	for enum, residue in enumerate(bchla_residues):
		#xyz = "["
		
		#xyz_file = "140\n"
		#xyz_file += f"Bchla {enum+1}\n"
		
	
		#for atom in residue.atoms():
		#	xyz += f"[{atom.element.symbol}, {str(positions[atom.index][0])[:-2]}, {str(positions[atom.index][1])[:-2]}, {str(positions[atom.index][2])[:-2]}], "
		#	xyz_file += f"{atom.element.symbol} {str(positions[atom.index][0])[:-2]} {str(positions[atom.index][1])[:-2]} {str(positions[atom.index][2])[:-2]}\n"
		
		#xyz = xyz[:-2]
		#xyz += "]"

		#print(xyz_file, file=open(f"Bchla_{enum+1}.xyz", 'w'))
		#result = run_qcore(xyz)

		result = read_monomer_result(enum)
		bchla_results.append(result)


	for enum1, residue1 in enumerate(bchla_residues):
		for enum2, residue2 in enumerate(bchla_residues):
			if enum1 >= enum2:
				continue

			#print(f"BChla {enum1+1} -- BChla {enum2+1}")
			#dimer = "["
			
			#xyz_file = "280\n"
			#xyz_file += f"Bchla {enum1+1} -- Bchla {enum2+1}\n"

			#for atom in residue1.atoms():
			#	dimer += f"[{atom.element.symbol}, {str(positions[atom.index][0])[:-2]}, {str(positions[atom.index][1])[:-2]}, {str(positions[atom.index][2])[:-2]}], "
			#	xyz_file += f"{atom.element.symbol} {str(positions[atom.index][0])[:-2]} {str(positions[atom.index][1])[:-2]} {str(positions[atom.index][2])[:-2]}\n"
			
			#for atom in residue2.atoms():
			#	dimer += f"[{atom.element.symbol}, {str(positions[atom.index][0])[:-2]}, {str(positions[atom.index][1])[:-2]}, {str(positions[atom.index][2])[:-2]}], "
			#	xyz_file += f"{atom.element.symbol} {str(positions[atom.index][0])[:-2]} {str(positions[atom.index][1])[:-2]} {str(positions[atom.index][2])[:-2]}\n"
			
			#print(xyz_file, file=open(f"Bchla_{enum1+1}_Bchla_{enum2+1}.xyz", 'w'))
			#xyz_dimer = dimer[:-2]
			#xyz_dimer += "]"

			#dimer_result = run_qcore(xyz_dimer)
			
			distance = get_distance(residue1, residue2, positions)
			dimer_result = read_dimer_result(enum1, enum2)
			charges = frenkel_hamiltonian_energy(residue1, residue2, positions, bchla_results[enum1], bchla_results[enum2], "charges")
			raw = frenkel_hamiltonian_energy(residue1, residue2, positions, bchla_results[enum1], bchla_results[enum2], "raw_dipole")
			orient = frenkel_hamiltonian_energy(residue1, residue2, positions, bchla_results[enum1], bchla_results[enum2], "NaNc_dipole")
			
			distances.append(distance)
			dimer_energies.append(dimer_result.energy)
			charge_coupling_energies.append(charges)
			dipole_coupling_energies.append(raw)
			orient_coupling_energies.append(orient)

			print(enum1, enum2)

	distances = np.array(distances)
	dimer_energies = np.array(dimer_energies)
	charge_coupling_energies = np.array(charge_coupling_energies)
	dipole_coupling_energies = np.array(dipole_coupling_energies)
	orient_coupling_energies = np.array(orient_coupling_energies)

	distances.dump("distances.pkl")
	dimer_energies.dump("dimer_energies.pkl")
	charge_coupling_energies.dump("charge_coupling_energies.pkl")
	dipole_coupling_energies.dump("dipole_coupling_energies.pkl")
	orient_coupling_energies.dump("orient_coupling_energies.pkl")

