import numpy as np

import data_objects

def angle(vec1, vec2):
    num = np.dot(vec1, vec2)
    dom = np.linalg.norm(vec1) * np.linalg.norm(vec2)
    
    angle = np.rad2deg(np.arccos(num/dom))
    
    return min(angle, 180-angle)


def dipole_coupling(tdm1, tdm2, structure1, structure2):
    r_vec = structure1.center - structure2.center
    r_mag = np.linalg.norm(r_vec)
    
    Nac1 = structure1.Na_Nc * (np.linalg.norm(tdm1)/np.linalg.norm(structure1.Na_Nc))
    Nac2 = structure2.Na_Nc * (np.linalg.norm(tdm2)/np.linalg.norm(structure2.Na_Nc))
        
    v = (np.dot(tdm1, tdm2)/(r_mag ** 3)) - 3*((np.dot(tdm1, r_vec) * np.dot(tdm2, r_vec))/(r_mag **5))
    
    return v

def charge_coupling(tcs1, tcs2, bchla1_structure, bchla2_structure):
    result = 0
    
    for n_i, i in enumerate(bchla1_structure.coords):
        for n_j, j in enumerate(bchla2_structure.coords):
            r = np.linalg.norm(i-j)
            
            result += tcs1[n_i] * tcs2[n_j] / r
            
    return float(result)

def switch_coupling_func(coupling):
    if coupling == "dipole":
        return lambda tdm1, tdm2, s1, s2 : dipole_coupling(tdm1, tdm2, s1, s2)
    
    elif coupling == "charges":
        return lambda tcs1, tcs2, s1, s2 : charge_coupling(tcs1, tcs2, s1, s2)

def state_hamiltonian(results, structures, coupling):
    assert(len(results) == len(structures))
    
    assert(coupling == "dipole" or coupling == "charges")
        
    n_sites = len(structures)
    n_states = n_sites + 1
    
    hamiltonian = np.zeros((n_states, n_states))
        
    #i,j is in the Hamiltonian basis - 0 is ground state, 1,... are excited states
    #m,n is in site basis
    
    for i in range(n_states):
        for j in range(n_states):
            if j < i:
                continue
        
            if i == j:
                #energy element
                if i == 0 and j == 0:
                    #ground energy with ground n-th chromophore
                    energy = sum([x.total_energy for x in results])

                    interaction = 0
                    for m in range(n_sites):
                        for n in range(n_sites):
                            if m >= n:
                                continue
                            
                            else:
                                if coupling == "dipole":
                                    interaction += dipole_coupling(results[m].molecular_dipole, results[n].molecular_dipole, structures[m], structures[n])
                                elif coupling == "charges":
                                    interaction += charge_coupling(results[m].partial_charges, results[n].partial_charges, structures[m], structures[n])
                                    
                    hamiltonian[i][j] = energy + interaction

                else:
                    #excited m-th energy with ground n-th chromophore
                    energy = sum([x.total_energy for x in results]) + results[i-1].transition_energy
                    interaction = 0
                    for n in range(n_sites):
                        if n+1 == i:
                            continue
                        else:
                            if coupling == "dipole":
                                interaction += dipole_coupling(results[i-1].excited_molecular_dipole, results[n].molecular_dipole, structures[i-1], structures[n])
                            elif coupling == "charges":
                                interaction += charge_coupling(results[i-1].excited_partial_charges, results[n].partial_charges, structures[i-1], structures[n])
                            
                    for m in range(n_sites):
                        for n in range(n_sites):
                            if m >= n:
                                continue
                            elif m+1 == i or n+1 == i:
                                continue
                            else: 
                                if coupling == "dipole":
                                    interaction += dipole_coupling(results[m].molecular_dipole, results[n].molecular_dipole, structures[m], structures[n])
                                elif coupling == "charges":
                                    interaction += charge_coupling(results[m].partial_charges, results[n].partial_charges, structures[m], structures[n])
                    
                    hamiltonian[i][j] = energy + interaction

            elif j > i:
                #coupling element
                if i == 0:
                    #coupling of ground state to all excited states bar the i-th
                    for m in range(n_sites):
                        if m == i:
                            continue
                        else:
                            if coupling == "dipole":
                                hamiltonian[i][j] += dipole_coupling(results[i].molecular_dipole, results[m].transition_dipole, structures[i], structures[m])
                            elif coupling == "charges":
                                hamiltonian[i][j] += charge_coupling(results[i].partial_charges, results[m].transition_charges, structures[i], structures[m])
            
                else:
                    #coupling of excited i-th chromophore to j-th chromophore
                    if coupling == "dipole":
                        hamiltonian[i][j] += dipole_coupling(results[i-1].transition_dipole, results[j-1].transition_dipole, structures[i-1], structures[j-1])
                    elif coupling == "charges":
                        hamiltonian[i][j] += charge_coupling(results[i-1].transition_charges, results[j-1].transition_charges, structures[i-1], structures[j-1])
                    
    return np.tril(hamiltonian.T) + np.triu(hamiltonian, k=1)

def excitation_hamiltonian(results, structures, coupling):
    assert(len(results) == len(structures))
        
    n_sites = len(structures)
    
    hamiltonian = np.zeros(n_sites, n_sites)
    
    for i in range(n_sites):
        for j in range(n_sites):
            if j < i:
                continue
                
            if i == j:
                hamiltonian[i][j] = results[i].transition_energy
                
            else:
                hamiltonian[i][j] = dipole_coupling(results[i].transition_dipole, results[j].transition_dipole, structures[i], structures[j])
                
    return hamiltonian
                

def switch_coupling(results, structures, coupling, state):
    if state:
        return state_hamiltonian(results, structures, coupling)
    else:
        return excitation_hamiltonian(results, structures, coupling)
    
def run_Frenkel_Hamiltonian(results, structures, coupling, state):
    hamiltonian = switch_coupling(results, structures, coupling, state)
    tdms_angle = angle(results[0].transition_dipole, results[1].transition_dipole)

    eigvals, eigvecs = np.linalg.eig(hamiltonian)
    
    return data_objects.FrenkelResult(hamiltonian, tdms_angle, eigvals, eigvecs)