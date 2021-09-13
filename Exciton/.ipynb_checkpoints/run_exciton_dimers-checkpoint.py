import subprocess
import os
import json

QCORE_PATH = os.environ["QCORE_PATH"]

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

    json_results = json.loads(json_run.stdout)
   
    print(json_results)

    #    arr = np.frombuffer(bytes.fromhex(" ".join([format(n, "02x") for n in data["water"]["transition_dipole"]["data"]["bytes"]]), dtype=data["water"]["transition_dipole"]["dtype"])
                        
    #print(json_results)

    

def get_current_results():
    with open("exciton_data.json") as f:
        return json.load(f)


if __name__ == "__main__":
    #json_file = get_current_results()
    
    for frame in range(1, 49951, 100):
        for i in range(1, 28):
            for j in range(i+1, 28):

                monomerA=f"trunc_bchla_{i}_frame_{frame}.xyz"
                monomerB=f"trunc_bchla_{j}_frame_{frame}.xyz"
                
                qcore_string = f"res := excitons(structure(file = '../monomer_xyzs/{monomerA}') structure(file = '../monomer_xyzs/{monomerB}') use_xtb = true hamiltonian = states)"
                
                run_qcore(qcore_string)
                
                exit()
                
                
                
                
                
                