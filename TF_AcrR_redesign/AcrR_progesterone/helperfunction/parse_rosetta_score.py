import os
import json
import pandas as pd
from multiprocessing import Pool, cpu_count
from pathlib import Path
import eval_mut_comb
import pyrosetta

# Initialize PyRosetta
pyrosetta.init("-mute all -ignore_unrecognized_res -ex1 -ex2aro -use_input_sc")

def process_single_pdb(args):
    """
    Wrapper function for parallel processing with error handling
    
    Args:
        args (tuple): Contains the sub_path and prot_name.
    Returns:

        tuple: (name_prefix, rosetta_energy, error_message)
        name_prefix (str): The name prefix derived from the sub_path.
        rosetta_energy (float or None): The Rosetta energy score, or None if an error occurs.
        error_message (str or None): Error message if an error occurs, otherwise None.

        updating 20250617
    
    """
    sub_path, prot_name = args
    try:
        name_prefix = os.path.basename(sub_path)
        predict_path = os.path.join(sub_path, "predictions")
        affinity_files = os.listdir(predict_path)
        affinity_files_real = os.path.join(predict_path, affinity_files[0])
        
        pdb_path = os.path.join(affinity_files_real, f"{affinity_files[0]}_model_0_cleanH.pdb")
        if not os.path.exists(pdb_path):
            return (name_prefix, None, "No PDB file found")
        
        rosetta_energy = eval_mut_comb.get_interface_scores(pdb_path)
        print(f"Processed {name_prefix} with energy: {rosetta_energy}")
        return (name_prefix, rosetta_energy, None)
    
    except Exception as e:
        return (name_prefix, None, str(e))

def prase_rosetta_score(pdb_boltz_folder: str, prot_name: str):
    """
    Parse Rosetta scores with multiprocessing support
    Args:
        pdb_boltz_folder (str): The path to the folder containing Boltz predictions.
        prot_name (str): The name of the protein, e.g. '1flm'.
    Returns:
        None: Saves the results to a JSON file and a CSV file with errors.

        updating 20250617
    
    """
    subfolders = [f.path for f in os.scandir(pdb_boltz_folder) if f.is_dir()]
    score_all = {}
    error_df = pd.DataFrame(columns=["prot", "error"])
    
    # Prepare tasks for multiprocessing
    tasks = [(sub_path, prot_name) for sub_path in subfolders]
    
    # Use 75% of available CPUs to avoid overloading
    num_processes = max(1, int(cpu_count() * 0.5))
    
    with Pool(processes=num_processes) as pool:
        results = pool.map(process_single_pdb, tasks)
    
    # Process results
    for name_prefix, energy, error in results:
        if error:
            print(f"Error processing {name_prefix}: {error}")
            error_df = pd.concat([
                error_df,
                pd.DataFrame({"prot": [name_prefix], "error": [error]})
            ], ignore_index=True)
        else:
            score_all[name_prefix] = energy
    
    # Save outputs
    error_path = os.path.join(pdb_boltz_folder, f"{prot_name}_pdb_errors.csv")
    error_df.to_csv(error_path, index=False)
    
    json_path = os.path.join(pdb_boltz_folder, f"{prot_name}_rosetta_energy.json")
    with open(json_path, 'w') as f:
        json.dump(score_all, f, indent=4)
    
    print(f"Processed {len(score_all)}/{len(subfolders)} structures successfully")

if __name__ == "__main__":
    pdb_boltz_folder = "/home/hdwang/protein_sensor_hd/1flm/1.Parse_result/1flm_interface_pdb_boltz_apo"
    prot_name = "1flm"
    prase_rosetta_score(pdb_boltz_folder, prot_name)