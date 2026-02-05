__version__ = "0.1.7"
# updating by hdwang on 20250715
# formatting boltz yaml 20260205
import os
import sys


# Allow running both:
#  - as a module:  python -m helperfunction.boltzutils
#  - as a script:  python /path/.../helperfunction/boltzutils.py
if __package__ in (None, ""):
    # Add the parent directory (which contains "helperfunction") to sys.path
    pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if pkg_root not in sys.path:
        sys.path.insert(0, pkg_root)
# Import eval_mut_comb depending on how we're invoked
if __package__:
    from . import eval_mut_comb
else:
    from helperfunction import eval_mut_comb

from helperfunction import eval_mut_comb
import yaml
from Bio import SeqIO
import json
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from pathlib import Path




# Function to convert FASTA to Boltz2 YAML with exact template matching
# Updating by hdwang on 20250701


###################################################################################################################
##############################section1 for Boltz PDB preparation and alignment####################################
##################################################################################################################

def fasta_to_boltz_yaml(fasta_file: str, output_yaml: str, ligand_smiles: Optional[str] = None, binder_id: Optional[str] = None) -> None:
    """
    Strictly convert FASTA to Boltz2 YAML with exact template matching.
    
    Args:
        fasta_file (str): Input FASTA file path.
        output_yaml (str): Output YAML file path.
        ligand_smiles (str, optional): SMILES string of the ligand.
        binder_id (str, optional): ID of the protein chain that binds the ligand (e.g., "C").
    """

    ''' 

        From a fasta file, we read the sequences, generate boltz header

        fasta format:
        >protein|name=chain_A
        MLPGTFFEVLKNEGVVAIATQGEDGPHLVNTWNSYLKVLDGNRIVVPVGGYHKTEANVARDERVLMTLGSRKVAGRNGPGTGFLIHGSAAFRTDGPEFEAIARFKWARAALVITVVSAEQTL
        >protein|name=chain_B
        MLPGTFFEVLKNEGVVAIATQGEDGPHLVNTWNSYLKVLDGNRIVVPVGGYHKTEANVARDERVLMTLGSRKVAGRNGPGTGFLIHGSAAFRTDGPEFEAIARFKWARAALVITVVSAEQTL



        e.g. Format:
        version: 1  # Optional, defaults to 1
        sequences:
        - protein:
            id: A
            sequence: MLPGTFFEVLKNEGVVAIATQGEDGPHLVNTWNSYLKVLDGNRIVVPVGGYHKTEANVARDERVLMTLGSRKVAGRNGPGTGFLIKGSAAFRTDGPEFEAIARFKWARAALVITVVSAEQTL
        - protein:
            id: B
            sequence: MLPGTFFEVLKNEGVVAIATQGEDGPHLVNTWNSYLKVLDGNRIVVPVGGYHKTEANVARDERVLMTLGSRKVAGRNGPGTGFLIKGSAAFRTDGPEFEAIARFKWARAALVITVVSAEQTL

        - ligand:
            id: C
            smiles: C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C

        properties:
        - affinity:
            binder: C

    '''

    # Initialize data structure
    data = {
        "version": 1,
        "sequences": [],
        "properties": [] if (ligand_smiles and binder_id) else None # switch to decide if predict affinity
    }
    
    # Parse FASTA and add protein sequences
    for record in SeqIO.parse(fasta_file, "fasta"):
        chain_id = record.id.split("=")[-1].strip()[-1]  # Extract "A" from "chain_A"
        data["sequences"].append({
            "protein": {
                "id": chain_id,
                "sequence": str(record.seq)
            }
        })
    
    # Add ligand if provided
    if ligand_smiles:
        data["sequences"].append({
            "ligand": {
                "id": "C",  # Hardcoded to match template
                "smiles": ligand_smiles
            }
        })
        
        # Add affinity property
        if binder_id == "C":  # Enforce binder_id="C" to match template
            data["properties"] = [{
                "affinity": {
                    "binder": "C"
                }
            }]
    
    # Custom YAML dumper for perfect formatting
    class IndentDumper(yaml.SafeDumper):
        def increase_indent(self, flow=False, indentless=False):
            return super().increase_indent(flow, False)  # Force indent after sequences/properties
    
    # Write YAML with exact indentation
    with open(f"yaml4boltz/{output_yaml}", 'w') as f:
        # create the directory if it does not exist
        os.makedirs(os.path.dirname(f.name), exist_ok=True)
        yaml.dump(
            data,
            f,
            Dumper=IndentDumper,
            sort_keys=False,
            indent=2,
            default_flow_style=False,
            width=1000
        )
        f.write("\n")  # Final newline to match template

# Example usage (generates identical output to your template)


def fasta_to_boltz_yaml_affinity_off(
    fasta_file,
    output_yaml,
    folder,
    ligand_smiles=None,
    binder_id=None,
    ligand_count=1,         # <— control how many ligands to add (0, 1, 2, …)
    ligand_ids=None         # <— optional explicit IDs, e.g., ["C"] or ["C","D"]
):
    """
    Convert FASTA to Boltz2 YAML. If ligand_smiles is provided, add ligand_count ligands.
    If ligand_ids is None, IDs are auto-generated starting from 'C' (skipping used protein IDs).
    """
    os.makedirs(folder, exist_ok=True)

    data = {
        "version": 1,
        "sequences": [],
        # "properties": []  # keep commented to keep affinity off
    }

    # Proteins
    for record in SeqIO.parse(fasta_file, "fasta"):
        chain_id = record.id.split("=")[-1].strip()[-1]  # adjust if your FASTA headers differ
        data["sequences"].append({
            "protein": {
                "id": chain_id if not binder_id else binder_id,  # use binder_id if you want to override
                "sequence": str(record.seq)
            }
        })

    # Ligands
    if ligand_smiles:
        used_ids = {e["protein"]["id"] for e in data["sequences"] if "protein" in e}

        # Normalize ligand_smiles to a list
        if isinstance(ligand_smiles, str):
            smiles_list = [ligand_smiles] * max(1, int(ligand_count))
        else:
            smiles_list = list(ligand_smiles)
            ligand_count = len(smiles_list)

        # Generate/check ligand_ids
        if ligand_ids is None:
            ligand_ids = []
            code = ord("C")
            while len(ligand_ids) < ligand_count:
                candidate = chr(code)
                code += 1
                if candidate in used_ids:
                    continue
                ligand_ids.append(candidate)
        else:
            if len(ligand_ids) != ligand_count:
                raise ValueError(f"ligand_ids length ({len(ligand_ids)}) must equal ligand_count ({ligand_count})")

        # Append ligands
        for lig_id, smi in zip(ligand_ids, smiles_list):
            data["sequences"].append({
                "ligand": {
                    "id": lig_id,
                    "smiles": smi
                }
            })
    
    # add pocket constraint for cofolding
    '''constraints:
    - pocket:
        binder: B1
        contacts: [ [ A1, 829 ], [ A1, 138 ] ]'''
    ############### add pocket constraint section ##################




    # Write YAML
    out_path = os.path.join(folder, output_yaml)
    with open(out_path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)
    return out_path



def prase_boltz_affinity_results(pdb_boltz_folder: str, prot_name: str):
    """
    Parse Boltz affinity results from a given folder and save them to a CSV file.
    
    Args:
        pdb_boltz_folder (str): The path to the folder containing Boltz predictions.
        prot_name (str): The name of the protein, e.g. '1flm'.
        
    Returns:
        writes a CSV file with the results in the format:
        prot_name, affinity_pred_value, affinity_probability_binary, affinity_pred_value1, affinity_pred_value2, affinity_probability_binary1, affinity_probability_binary2
    """
    # get subfolders in the folder
    subfolders = [f.path for f in os.scandir(pdb_boltz_folder) if f.is_dir()]

    
    # direct to the predictions folder
    error_df = pd.DataFrame(columns=["prot", "error"])
    score_all = {}  # to collect all scores
    for sub_path in subfolders:
    
        name_prefix = sub_path.split("/")[-1]  # get the name of the subfolder, e.g. 1flm_hcy12_87
    
    # crate a dictionary to store the scores for this protein
        if name_prefix not in score_all:
            score_all[name_prefix] = {}

        # read affinity json file which starts with "affinity_"
        predict_path = os.path.join(sub_path, "predictions")
        print(predict_path)
        affinity_files = os.listdir(predict_path)


        try:
            affinity_files_real = os.path.join(predict_path, affinity_files[0])
            
            print(affinity_files_real)

            affinity_files = [f for f in os.listdir(affinity_files_real) if f.startswith("affinity_") and f.endswith(".json")]
            print(affinity_files)
            # get the absolute path of the first affinity file
            if not affinity_files:
                print(f"No affinity files found in {affinity_files_real}. Skipping...")
                continue
            affinity_file_path = os.path.join(affinity_files_real, affinity_files[0])

            


            # read the json data
        except IndexError:
            print(f"No affinity files found in {predict_path}. Skipping...")
            error_df = error_df.append({"prot": name_prefix, "error": "No affinity files found"}, ignore_index=True)
            error_df.to_csv(os.path.join(pdb_boltz_folder, f"{prot_name}_affinity_errors.csv"), index=False)
            continue

            
            # Read the JSON data
        with open(affinity_file_path, 'r') as f:

            data = json.load(f)
            affinity_value = data["affinity_pred_value"]
            affinity_probability_binary = data["affinity_probability_binary"]
            affinity_pred_value1 = data["affinity_pred_value1"]
            affinity_pred_value2 = data["affinity_pred_value2"]
            affinity_probability_binary1 = data["affinity_probability_binary1"]
            affinity_probability_binary2 = data["affinity_probability_binary2"]

            #collect all values if exists 
            
            score_all[name_prefix]["affinity_pred_value"] = affinity_value
            score_all[name_prefix]["affinity_probability_binary"] = affinity_probability_binary
            score_all[name_prefix]["affinity_pred_value1"] = affinity_pred_value1
            score_all[name_prefix]["affinity_pred_value2"] = affinity_pred_value2
            score_all[name_prefix]["affinity_probability_binary1"] = affinity_probability_binary1  
            score_all[name_prefix]["affinity_probability_binary2"] = affinity_probability_binary2

    # generate a csv file with the results
    
    csv_file = os.path.join(pdb_boltz_folder,f"{prot_name}_affinity_results.csv")

    df_scores = pd.DataFrame.from_dict(score_all, orient='index')
    df_scores.to_csv(csv_file, index_label='prot_name')
    # or write them into json file
    json_file = os.path.join(pdb_boltz_folder, f"{prot_name}_affinity_results.json")
    with open(json_file, 'w') as f:
        json.dump(score_all, f, indent=4)

CONF_KEYS = [
    "confidence_score", "ptm", "iptm", "ligand_iptm", "protein_iptm",
    "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde",
]


def _safe_mean(vals: List[float]) -> Optional[float]:
    if not vals:
        return None
    arr = np.asarray(vals, dtype=float)
    if arr.size == 0:
        return None
    return float(np.nanmean(arr))

def prase_boltz_confidence_results(pdb_boltz_folder: str, prot_name: str, wt_pdb: str):
    """
    Parse confidence_* results under each subfolder's 'predictions' directory.
    Averages metrics across all matching confidence_* entries per subfolder.

    Outputs:
      - {prot_name}_confidence_results.csv (aggregated means per subfolder)
      - {prot_name}_confidence_results.json (same content as a dict)
      - {prot_name}_confidence_errors.csv (rows that failed/empty)
    """
    pdb_boltz_folder = str(pdb_boltz_folder)
    subfolders = [f.path for f in os.scandir(pdb_boltz_folder) if f.is_dir()]
    subfolders = sorted(subfolders)

    error_df = pd.DataFrame(columns=["prot", "error"])
    score_all: Dict[str, Dict] = {}

    for sub_path in subfolders:
        name_prefix = Path(sub_path).name  # e.g., 1flm_hcy12_87
        score_all.setdefault(name_prefix, {})
        # under sub_path, there is dir using the name,
        predictions_path_upper = os.path.join(sub_path, "predictions")
        working_path = os.listdir(predictions_path_upper)
        predict_path = Path(predictions_path_upper) / working_path[0]

        print(f"Looking for predictions in: {predict_path}")

        if not os.path.isdir(predictions_path_upper):
            msg = f"Missing predictions/ in {predictions_path_upper}"
            print(f"[WARN] {msg}")
            error_df = pd.concat([error_df, pd.DataFrame({"prot": [name_prefix], "error": [msg]})],
                                 ignore_index=True)
            continue

        # Find all entries starting with "confidence_" inside predictions/
        entries = [e for e in os.scandir(predict_path) if e.name.startswith("confidence_")]
        if not entries:
            msg = f"No entries starting with 'confidence_' in {predict_path}"
            print(f"[WARN] {msg}")
            error_df = pd.concat([error_df, pd.DataFrame({"prot": [name_prefix], "error": [msg]})],
                                 ignore_index=True)
            continue

        print(f"processing {predict_path}, found {len(entries)} 'confidence_' entries")

        # Accumulate metrics across all model JSONs
        accum = {k: [] for k in CONF_KEYS}
        used_models = 0

        for e in sorted(entries, key=lambda x: x.name):
            json_path = None

            if e.is_file() and e.name.endswith(".json"):
                # e.g., predictions/confidence_..._model_0.json
                json_path = e.path
            elif e.is_dir():
                # Prefer a JSON starting with confidence_ inside this directory
                files = [f for f in os.listdir(e.path) if f.endswith(".json")]
                prefer = [f for f in files if f.startswith("confidence_")]
                pick = prefer[0] if prefer else (files[0] if files else None)
                if pick:
                    json_path = os.path.join(e.path, pick)

            if not json_path:
                print(f"[WARN] No JSON found for {e.name}, skipping.")
                continue

            try:
                with open(json_path, "r") as f:
                    data = json.load(f)
                for k in CONF_KEYS:
                    if k in data:
                        try:
                            accum[k].append(float(data[k]))
                        except Exception:
                            pass
                used_models += 1
            except Exception as ex:
                print(f"[WARN] Failed reading {json_path}: {ex}")

        if used_models == 0:
            msg = f"No readable confidence JSONs under entries in {predict_path}"
            print(f"[WARN] {msg}")
            error_df = pd.concat([error_df, pd.DataFrame({"prot": [name_prefix], "error": [msg]})],
                                 ignore_index=True)
            continue

        avg_metrics = {k: _safe_mean(v) for k, v in accum.items()}
        avg_metrics["n_models"] = int(used_models)

        score_all[name_prefix].update(avg_metrics)
        print(f"Averages for {name_prefix} over {used_models} model(s): "
              f"confidence_score={avg_metrics.get('confidence_score')}")

    # Write aggregated results
    out_csv = os.path.join(pdb_boltz_folder, f"{prot_name}_confidence_results.csv")
    df_scores = pd.DataFrame.from_dict(score_all, orient="index")
    df_scores.index.name = "prot_name"
    df_scores.to_csv(out_csv)

    out_json = os.path.join(pdb_boltz_folder, f"{prot_name}_confidence_results.json")
    with open(out_json, "w") as f:
        json.dump(score_all, f, indent=4)

    out_err = os.path.join(pdb_boltz_folder, f"{prot_name}_confidence_errors.csv")
    if not error_df.empty:
        error_df.to_csv(out_err, index=False)
    else:
        # ensure old error file doesn't mislead
        if os.path.exists(out_err):
            try:
                os.remove(out_err)
            except OSError:
                pass

    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_json}")
    if not error_df.empty:
        print(f"Wrote: {out_err}")
    else:
        print("No errors encountered.")



###################################################################################################################
##############################section2 for Boltz PDB preparation and alignment####################################
##################################################################################################################



import os
import pandas as pd
from pymol import cmd,finish_launching



def clean_boltzpdb(pdb_path: str):
        
    """
    Prepare a Boltz PDB by removing ligands and adding hydrogens.
    Outputs (if present):
    - <basename>_ligandC.pdb : ligand with resn LIG and chain C
    - <basename>_ligandD.pdb : ligand with resn LIG and chain D
    - <basename>_cleanH.pdb  : molecule without ligands, with hydrogens added

    Args:
        pdb_path (str): Path to the input PDB file.

    Returns:
        dict: Paths to created files keyed by 'ligandC', 'ligandD', 'cleanH' (only keys that were created).
    """
    finish_launching(['pymol', '-cq'])  # headless + quiet

    base, _ = os.path.splitext(pdb_path)
    out_clean = f"{base}_cleanH.pdb"
    out_ligC = f"{base}_ligandC.pdb"
    out_ligD = f"{base}_ligandD.pdb"

    created = {}

    try:
        cmd.reinitialize()  # start clean
        cmd.load(pdb_path, 'mol')

        # Save ligand C if present
        nC = cmd.count_atoms("resn LIG and chain C")
        if nC > 0:
            cmd.save(out_ligC, "resn LIG and chain C")
            created['ligandC'] = out_ligC
        else:
            print("No ligand C found; skipping ligandC output.")

        # Save ligand D only if present
        nD = cmd.count_atoms("resn LIG and chain D")
        if nD > 0:
            cmd.save(out_ligD, "resn LIG and chain D")
            created['ligandD'] = out_ligD
        else:
            print("No ligand D found; skipping ligandD output.")

        # Remove all ligands, add hydrogens to the remaining structure, and save
        cmd.remove("resn LIG")
        cmd.h_add("mol")  # add hydrogens to the cleaned molecule
        cmd.save(out_clean, "mol")
        created['cleanH'] = out_clean

    finally:
        cmd.delete("all")

def parse_boltzpdb_rmsd(pdb_path: str,refpdb_path: str,output_folder: str = "aligned_pdb_boltz"):
    """
    Prepare a Boltz PDB file by removing ligands and adding hydrogens.
    
    Args:
        pdb_path (str): Path to the input PDB file.
        refpdb_path (str): Path to the reference PDB file for alignment.
        output_folder (str): Folder to save the aligned PDB file.
    Returns:
        alignment_df (pd.DataFrame): DataFrame containing the protein name and RMSD score.
    """
    
    
    finish_launching(['pymol', '-c'])  # '-c' ensures it runs without GUI
    print(f"Loading {pdb_path} and {refpdb_path} for alignment...")
    cmd.delete("molecule")  # Clear previous molecule if exists
    cmd.load(pdb_path, 'molecule')
    cmd.load(refpdb_path, 'reference')
    
    alignment = cmd.align("molecule and name CA", "reference and name CA", cycles=10)
    print(f"Aligned {pdb_path} to reference, ")
    print(f"Alignment score_{pdb_path}:", alignment[0])


    #cmd.h_add('molecule')  # Add hydrogens

    # revise the output file name, add  alignment rmsd behind the name
    rmsd = alignment[0]
    name = os.path.basename(pdb_path).split('.')[0]  # Get the base name without extension
    output_folder = os.path.join(output_folder, name)
    os.makedirs(output_folder, exist_ok=True)  # Create output folder if it doesn't exist
    output_file = os.path.join(output_folder, f'{name}_rmsd_{rmsd:.2f}.pdb')
    cmd.save(output_file, "molecule")

    # record all alignment score in a dataframe
    alignment_df = pd.DataFrame({
        'prot_name': [name],
        'rmsd': [rmsd]
    })
    
    return alignment_df

def prase_boltz_pdb_results(pdb_boltz_folder: str, prot_name: str, ref_pdb: str):

    """
    Parse Boltz pdb from given folder: commly I put them all in named: pdbid_pdb_boltz 
    
    Args:
        pdb_boltz_folder (str): The path to the folder containing Boltz predictions.
        prot_name (str): The name of the protein, e.g. '1flm'.
        
    Returns:
    """

    # get subfolders in the folder
    subfolders = [f.path for f in os.scandir(pdb_boltz_folder) if f.is_dir()]
    
    # direct to the predictions folder
    error_df = pd.DataFrame(columns=["prot", "error"])
    rmsd_df = pd.DataFrame(columns=["prot_name", "rmsd"])  # to collect all rmsd scores
    for sub_path in subfolders:
    
        name_prefix = sub_path.split("/")[-1]  # get the name of the subfolder, e.g. 1flm_hcy12_87
    


        # read affinity json file which starts with "affinity_"
        predict_path = os.path.join(sub_path, "predictions")
        print(f"Processing predictions in {predict_path}")

        affinity_files = os.listdir(predict_path)
        # make sure there is at least one affinity file
        if not affinity_files:
            print(f"No affinity files found in {predict_path}. Skipping...")
            # concat error_df with the name_prefix and error message
            df_temp = pd.DataFrame({"prot": [name_prefix], "error": ["No affinity files found"]})
            error_df = pd.concat([error_df, df_temp], ignore_index=True)
            continue
        affinity_files_real = os.path.join(predict_path, affinity_files[0]) # it should output the folder path containing the pdb files
        
        pdb_path = os.path.join(affinity_files_real, f"{affinity_files[0]}_model_0.pdb")

        # parse the pdb file, remove ligands, add hydrogens, save as a new pdb file
        aliged_df = parse_boltzpdb_rmsd(
            pdb_path=pdb_path,
            refpdb_path=ref_pdb
        )
        
        # collect all alignment scores
        rmsd_df = pd.concat([rmsd_df, aliged_df], ignore_index=True)
    # generate a csv file with the results
    csv_file = os.path.join(pdb_boltz_folder, f"{prot_name}_rmsd_results.csv")
    rmsd_df.to_csv(csv_file, index=False)

    error_csv = os.path.join(pdb_boltz_folder, f"{prot_name}_error_results.csv")
    error_df.to_csv(error_csv, index=False)

def do_boltzpdb_clean(pdb_boltz_folder: str, prot_name: str):

    """
    Parse Boltz pdb from given folder: commly I put them all and named: pdbid_pdb_boltz 
    
    Args:
        pdb_boltz_folder (str): The path to the folder containing Boltz predictions.
        prot_name (str): The name of the protein, e.g. '1flm'.
        
    Returns:
    """

    # get subfolders in the folder
    subfolders = [f.path for f in os.scandir(pdb_boltz_folder) if f.is_dir()]
    
    # direct to the predictions folder
    error_df = pd.DataFrame(columns=["prot", "error"])
    rmsd_df = pd.DataFrame(columns=["prot_name", "rmsd"])  # to collect all rmsd scores
    for sub_path in subfolders:
    
        name_prefix = sub_path.split("/")[-1]  # get the name of the subfolder, e.g. 1flm_hcy12_87
    


        # read affinity json file which starts with "affinity_"
        try:
            predict_path = os.path.join(sub_path, "predictions")
        
            affinity_files = os.listdir(predict_path)

            affinity_files_real = os.path.join(predict_path, affinity_files[0]) # it should output the folder path containing the pdb files
            
            pdb_path = os.path.join(affinity_files_real, f"{affinity_files[0]}_model_0.pdb")
        
            # prepare the pdb file, remove ligands, add hydrogens, save as a new pdb file
            clean_boltzpdb(pdb_path)
        except Exception as e:
            print(f"Error processing {sub_path}: {e}")
            # concat error_df with the name_prefix and error message
            df_temp = pd.DataFrame({"prot": [name_prefix], "error": [str(e)]})
            error_df = pd.concat([error_df, df_temp], ignore_index=True)
            continue

    error_csv = os.path.join(pdb_boltz_folder, f"{prot_name}_error_results.csv")
    error_df.to_csv(error_csv, index=False)

        
# Example usage:
#input_folder = "/mnt/d/binder_finder_backup/plip_prase/25.06/test4boltz2/3cyl_yaml"

#run_boltz_batch_prediction_sequentially(input_folder)



import argparse

# -------- argument validators --------
def existing_dir(path: str) -> str:
    if os.path.isdir(path):
        return path
    raise argparse.ArgumentTypeError(f"directory not found: {path}")

def existing_file(path: str) -> str:
    if os.path.isfile(path):
        return path
    raise argparse.ArgumentTypeError(f"file not found: {path}")

# -------- command handlers --------
def cmd_rmsd(args: argparse.Namespace) -> None:
    prase_boltz_pdb_results(
        pdb_boltz_folder=args.pdb_boltz_folder,
        prot_name=args.prot_name,
        ref_pdb=args.ref_pdb,
    )
    print("RMSD parsing complete.")

def cmd_clean(args: argparse.Namespace) -> None:
    do_boltzpdb_clean(
        pdb_boltz_folder=args.pdb_boltz_folder,
        prot_name=args.prot_name,
    )
    print("Cleaning complete.")

def cmd_affinity(args: argparse.Namespace) -> None:
    prase_boltz_affinity_results(
        pdb_boltz_folder=args.pdb_boltz_folder,
        prot_name=args.prot_name,
    )
    print("Affinity parsing complete.")

def cmd_confidence(args: argparse.Namespace) -> None:
    prase_boltz_confidence_results(
        pdb_boltz_folder=args.pdb_boltz_folder,
        prot_name=args.prot_name,
        wt_pdb=args.wt_pdb,
    )
    print("Confidence parsing complete.")

# -------- parser builder --------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="boltz-utils",
        description="Utilities to parse/clean Boltz PDB predictions.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # rmsd
    p_rmsd = subparsers.add_parser(
        "rmsd",
        help="Parse Boltz PDBs and compute RMSD vs a reference PDB.",
    )
    p_rmsd.add_argument(
        "--pdb-boltz-folder",
        required=True,
        type=existing_dir,
        help="Path to the folder containing Boltz prediction subfolders.",
    )
    p_rmsd.add_argument(
        "--prot-name",
        required=True,
        help="Protein name tag used for output CSV filenames.",
    )
    p_rmsd.add_argument(
        "--ref-pdb",
        required=True,
        type=existing_file,
        help="Reference PDB file path for RMSD alignment.",
    )
    p_rmsd.set_defaults(func=cmd_rmsd)

    # clean
    p_clean = subparsers.add_parser(
        "clean",
        help="Clean Boltz PDBs (remove ligands, add hydrogens) across prediction folders.",
    )
    p_clean.add_argument(
        "--pdb-boltz-folder",
        required=True,
        type=existing_dir,
        help="Path to the folder containing Boltz prediction subfolders.",
    )
    p_clean.add_argument(
        "--prot-name",
        required=True,
        help="Protein name tag used for error CSV filename.",
    )
    p_clean.set_defaults(func=cmd_clean)

    # affinity
    p_affi = subparsers.add_parser(
        "affinity",
        help="Parse Boltz affinity results across prediction folders.",
    )
    p_affi.add_argument(
        "--pdb-boltz-folder",
        required=True,
        type=existing_dir,
        help="Path to the folder containing Boltz prediction subfolders.",
    )
    p_affi.add_argument(
        "--prot-name",
        required=True,
        help="Protein name tag used for output CSV filenames.",
    )
    p_affi.set_defaults(func=cmd_affinity)

    # confidence
    p_conf = subparsers.add_parser(
        "confidence",
        help="Parse Boltz confidence results across prediction folders.",
    )
    p_conf.add_argument(
        "--pdb-boltz-folder",
        required=True,
        type=existing_dir,
        help="Path to the folder containing Boltz prediction subfolders.",
    )
    p_conf.add_argument(
        "--prot-name",
        required=True,
        help="Protein name tag used for output CSV filenames.",
    )
    p_conf.add_argument(
        "--wt-pdb",
        required=True,
        type=existing_file,
        help="Wild-type PDB file path for pLDDT comparison.",
    )
    p_conf.set_defaults(func=cmd_confidence)

    return parser

# -------- entry point --------
def main(argv=None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()