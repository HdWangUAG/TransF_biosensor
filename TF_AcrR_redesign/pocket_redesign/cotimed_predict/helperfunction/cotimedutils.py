__version__ = "0.1.7"
# Updated by H.D.Wang, 2025.08.13
# Formatted 27.Jan.2026
# This script contains utility functions for inverse protein design using TIMED.

import os
import json
from pathlib import Path
from itertools import product
from collections import defaultdict
from typing import Dict, List

import pandas as pd
import numpy as np
import ampal
import ampal.geometry

import argparse

# ================================================================
#  Section 1: Protein structure utilities
# ================================================================



def residue_number_indexing(assembly: ampal.Assembly):
    
    """
    Assigns sequential residue IDs to an AMPAL assembly.
    """
    count = 0
    residue_number = []

    for chain in assembly:
        for residue in chain:
            count += 1
            residue.id = count
            residue_number.append(int(residue.id))


def generate_fasta(protein_path: str) -> Dict[int, str]:
    """
    Generates a residue index-to-AA dictionary from a PDB file.

    Parameters
    ----------
    protein_path : str
        Path to the PDB file.

    Returns
    -------
    dict
        {1: 'A', 2: 'C', 3: 'D', ...}
    """


    # Load the PDB file
    my_structure = ampal.load_pdb(protein_path)

    # get fasta sequence
        
    sequences = my_structure.sequences
    
    # extract the sequence from the list and convert it to a string
    fasta_seq = ''.join(sequences)
    
    # create a dictionary to store the residue number and the corresponding three letter code
    residue_dict = {}
    for i in range(0,len(fasta_seq)):
        residue_dict[i+1] = fasta_seq[i]
    
    return residue_dict

# for dimer get single chain
def get_single_chain(pdb_path: str, ) -> str:
    """
    Extracts the amino acid sequence from a PDB file and returns it as a string.
    Args:
        pdb_path (str): Path to the PDB file.
    Returns:
        str: Amino acid sequence extracted from the PDB file.
    """
    # Load the PDB file
    my_structure = ampal.load_pdb(pdb_path)

    # get fasta sequence
    sequences = my_structure.sequences
    
    # extract the sequence from the list and convert it to a string
    fasta_seq = sequences[0]
    
    return fasta_seq


# ================================================================
#  Section 2: Pocket residue extraction
# ================================================================


def get_pocket_residues(pocket_size: float) -> Dict[str, List[int]]:
    """
    Extract residues within a distance cutoff of ligands in all PDBs in cwd.

    Parameters
    ----------
    pocket_size : float
        Distance cutoff for pocket residues.

    Returns
    -------
    dict
        { "protein_chain": [residue indices, ...], ... }
    """
    ligands_of_interest = Path("ligand.csv").read_text().splitlines()
    files = [f for f in os.listdir() if f.endswith(".pdb")]
    backbone_atoms = {"C", "N", "O", "CA"}

    pocket_dict = defaultdict(list)

    for pdb_file in files:
        protein = ampal.load_pdb(pdb_file)
        residue_number_indexing(protein)

        for ligand in protein.get_ligands():
            if ligand.mol_code not in ligands_of_interest:
                continue

            chain_res_map = defaultdict(set)
            for atom_name in ligand.atoms:
                for atom in protein.get_atoms(ligands=False):
                    if atom.res_label not in backbone_atoms:
                        if ampal.geometry.distance(ligand[atom_name], atom) <= pocket_size:
                            chain_id = f"{Path(pdb_file).stem}_{atom.parent.parent.id}"
                            chain_res_map[chain_id].add(int(atom.parent.id))

            for chain, residues in chain_res_map.items():
                pocket_dict[chain].extend(residues)

    with open("pocket_residues.json", "w") as f:
        json.dump(pocket_dict, f, indent=4)

    return dict(pocket_dict)


# ================================================================
#  Section 3: Probability & label map processing
# ================================================================


def process_data(
    probabilities_file: str,
    label_map_file: str,
    output_file: str,
    pdb_name: str
) -> pd.DataFrame:
    """
    Processes and combines probabilities data with label map data,
    keeping only entries that match the given PDB name.

    Parameters
    ----------
    probabilities_file : str
        Path to the CSV file containing probabilities.
    label_map_file : str
        Path to the text file containing label mappings.
    output_file : str
        Desired name for the output CSV file.
    pdb_name : str
        PDB identifier to filter on (e.g., "3cyl_A").

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with probabilities.
    """
    alphabet = list("ACDEFGHIKLMNPQRSTVWY")

    # Read probabilities
    df_probs = pd.read_csv(probabilities_file, header=None, names=alphabet)

    # Read label map
    df_labels = pd.read_csv(
        label_map_file,
        header=None,
        names=["pdbinfo", "chain", "idx", "AAs"]
    )
    df_labels["pdbinfo_chain"] = df_labels["pdbinfo"] + "_" + df_labels["chain"]

    # Combine data
    df_combined = pd.concat([df_labels["pdbinfo_chain"], df_labels["idx"], df_probs], axis=1)

    # Filter only rows that contain the pdb_name
    df_filtered = df_combined[df_combined["pdbinfo_chain"].str.contains(pdb_name, case=False, na=False)]

    # Save filtered results
    df_filtered.to_csv(output_file, index=False)
    print(f"[INFO] Saved filtered data for '{pdb_name}' → {output_file}")

    return df_filtered

def get_amino_acid(
    pdb_chain: str, 
    residue_idx: int, 
    df: pd.DataFrame, 
    cutoff: float = 0.1
) -> List[str]:
    """
    Get probable amino acids for a residue above a probability cutoff.

    Parameters
    ----------
    pdb_chain : str
        Protein + chain identifier, e.g. "3cyl_A".
    residue_idx : int
        Residue index.
    df : pd.DataFrame
        DataFrame with probabilities.
    cutoff : float, optional
        Probability threshold. Defaults to 0.1.

    Returns
    -------
    list of str
        Amino acids whose predicted probability >= cutoff.
    """
    alphabet = list("ACDEFGHIKLMNPQRSTVWY")
    row = df[(df["pdbinfo_chain"] == pdb_chain) & (df["idx"] == residue_idx)]
    if row.empty:
        return []

    probs = row.iloc[0][alphabet]
    return [aa for aa in alphabet if probs[aa] >= cutoff]

# ================================================================
#  Section 4: Mutation design
# ================================================================


def generate_combination(prot_name: str, predicted_json_path: str) -> list:

    ''' 
        In cotimed, which mutate all residues the idx is from 1 to total number of residues,
        so we need to relabel the residue number, here I generate new dict to store all residues

        Arguments:
        prot_name: the name of the protein, e.g. '1flm_hcy12'
        predicted_json_path: the path to the predicted JSON file containing mutation predictions
        e.g. 

        residue_dict: a dictionary containing residue numbers and their corresponding amino acids
        e.g. {1: 'A', 2: 'C', 3: 'D', ...}
        This function generates all possible combinations of mutations for the specified protein.

        returns:
        permutation_matrix: a list of dictionaries, each dictionary represents a combination of mutations
        e.g. [{'1': 'A', '2': 'C', '3': 'D'}, {'1': 'A', '2': 'C', '3': 'E'}, ...]
    '''
    
    print(f'Processing {prot_name}')
    mutate_dict = {}
    
    with open(predicted_json_path, 'r') as f:
        all_predict = json.load(f)


    for key in all_predict:

        if prot_name in key:

            ''' no matter how many chains, we introduce all residues'''

            mutate_dict.update({k: v for k, v in all_predict[key].items() if v != []})
  
    # create permutation matrix to store the residue number that need to be mutated
    permutation_matrix = []

    # 1. sort the keys in the mutate_dict, and get the order of the residue position
    sorted_keys = sorted(mutate_dict.keys()) # create a list of position keys

    candidates_list = [mutate_dict[key] for key in sorted_keys]
 
    # 2.generate all possible combinations
    all_combinations = list(product(*candidates_list))


    # 3. add the residue number to the permutation matrix

    for i in all_combinations:
        temp = {}
        for j in range(len(i)):

            temp[sorted_keys[j]] = i[j]

        permutation_matrix.append(temp)

    return permutation_matrix


# mutate all residues, iterate the dict

def all_on_1chain(all_predict: List[Dict[str, str]], single_chain: str, new_name: str) -> Dict[str, str]:
    """
    Apply mutations across all residues of a single chain (for homodimers).
    """
    mutated_all = {}
    for i, mutation_set in enumerate(all_predict, start=1):
        mutated_seq = list(single_chain)
        for pos, aa in mutation_set.items():
            idx = int(pos) - 1
            if 0 <= idx < len(mutated_seq):
                mutated_seq[idx] = aa
        mutated_all[f"{new_name}_{i}"] = "".join(mutated_seq)
    return mutated_all

# ================================================================
#  Section 5: FASTA generation
# ================================================================


def write_fasta(seq_list,filename,save_path):
    """
    Write sequences to FASTA file with headers.
    """
    Path(save_path).mkdir(parents=True, exist_ok=True)
    filepath = Path(save_path) / filename
    
    with open(filepath, 'w') as f:
        header_prefix = '>protein|name=chain_'
        header_list = [header_prefix + chr(65 + i) for i in range(len(seq_list))]

        for seq,header in zip(seq_list,header_list):
            # add header
            f.write(header + '\n')
            # write sequence
            f.write(seq + '\n')


def write_fasta_from_dict(fasta_matrix, save_folder: str):
    '''
    
    input: fasta_matrix_folder: str, the path of the fasta sequence matrix
           save_folder: str, the path of the folder to save the fasta sequence
    output: save the fasta sequence to the save_folder
    description:
    1. read the fasta sequence matrix from result generated by all_on_1chain
    2. extract the key as filename
    3. extract the value as sequence
    4. save the sequence to a fasta file, as 2fnu is a dimer, so save the sequence twice
    5. check if the save folder exists if not create it
    6. write the fasta file
    
    revised by H.D.Wang 2025.04.12
    '''

        # open the fasta sequence matrix

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    
    for filename, sequence in fasta_matrix.items():
        # Save as <filename>.fasta
         
        write_fasta([sequence,sequence], filename,save_folder)  # Wrap sequence in a list for write_fasta
        
        

def do_unique_fasta(prot_name: str):
    '''
    This function is to generate a unique fasta file from the sequence library.
    Arguments:
        prot_name (str): The name of the protein, e.g. '1q12_hcy'.
    Returns:
        None
    '''
    current_path = os.getcwd()
    path_pdb = [f for f in os.listdir(current_path) if f.startswith(prot_name) and f.endswith('onechain_mut_all.json')]
    
    sequence_lib = {}
    existing_sequences = set()

    for path in path_pdb:
        with open(path, 'r') as f:
            data = json.load(f)

        sequence_lib.update(data)

    sequence_lib_unique = {k: v for k, v in sequence_lib.items() if v not in existing_sequences and not existing_sequences.add(v)}

    output_fasta_path = f'{prot_name}_unique_sequences'
    write_fasta_from_dict(sequence_lib_unique, output_fasta_path)
    print(f"Unique sequences saved to '{output_fasta_path}'.")
    print(f"Total unique sequences: {len(sequence_lib_unique)}")




def do_homodimer_pipeline(
    pdb_name: str,
    pocket_size: float,
    cutoff: float,
    probabilities_path: str = "cnn_timed_50percentmodel_cofactor-148-1.73.csv",
    label_map_path: str = "datasetmap.txt",
) -> None:
    '''
        The pipeline for sequence design, including:
        1. Process the generated CSV file from the timed design CSV, 'cnn_timed_50percentmodel_cofactor-148-1.73.csv'
        and 'datasetmap.txt', to get the predicted probabilities for each residue.
        2. Process ligand surrounding residues, 'pocket_residues.json',
        3. Filtering the residues based on the predicted probabilities and the cutoff value. 
        get a JSON file with the predicted sequence for each residue.e.g."1flm_hcy12_redesign_pocket_sequence.json"ArithmeticError
        4. Generate the mutation combinations based on the predicted JSON file.
        and here we assume the protein is a homodimer, so we need to mutate all residues in the single chain.
        Arguments:
            protein_name (str): The name of the protein, e.g. '1flm_hcy12'.
        
        Returns: generates a JSON file with the predicted sequence for each residue.


    '''

    # --- Section 1: Process generated CSV file from timed design CSV ---
    output_csv_name = f'{pdb_name}_predict_prob.csv'

    # Process data and get the DataFrame
    df_prob = process_data(probabilities_path, label_map_path, output_csv_name,pdb_name)
    
    # --- Section 2: Process the redesign sequence and JSON file from the selected sequence ---
    get_pocket_residues(pocket_size)
    current_directory = os.getcwd()

    # Load redesign sequence information from a JSON file

    json_input_path = os.path.join(current_directory, 'pocket_residues.json')
    with open(json_input_path, 'r') as f:
        redesign_sequence = json.load(f)

    predict_seq = {}
    for key_full_name in redesign_sequence.keys():

        predict_seq[key_full_name] = {}  # Initialize a sub-dictionary for each key_full_name
        for item_idx in redesign_sequence[key_full_name]:
            # Ensure item_idx is an integer for lookup
            predict_seq[key_full_name][item_idx] = get_amino_acid(key_full_name, int(item_idx), df_prob, cutoff)

    # Save the predicted sequence to a JSON file
    output_json_name = f'{pdb_name}_redesign_pocket_sequence.json'
    json_output_path = os.path.join(current_directory, output_json_name)

    with open(json_output_path, 'w') as f:
        json.dump(predict_seq, f, indent=4)
    print(f"Predicted sequence data saved to '{output_json_name}'.")

    # --- Section 3: Generate the mutatation combinations based on the predicted JSON file ---
    # Generate the mutation combinations
    process_pdb = pdb_name + '.pdb'
    redesign_sequence = generate_fasta(process_pdb)
    mutation_combinations = generate_combination(pdb_name, json_output_path)
    print(f"Generated {len(mutation_combinations)} mutation combinations.")
    print(f"Mutation combinations: {mutation_combinations}")
    
    # Save the mutation combinations to a CSV file
   
    single_chain = get_single_chain(process_pdb)
    print(f"Single chain for {pdb_name}: {single_chain}")

    mutated_dict = all_on_1chain(mutation_combinations, single_chain, pdb_name)
    

    # check if there has repeat sequences
    existing_sequences = set()
    sequence_lib_unique = {
        k: v for k, v in mutated_dict.items()
        if v not in existing_sequences and not existing_sequences.add(v)
    }

    # 4.1 Save unique sequences to JSON
    with open("mutated_homodimer_sequences.json", "w") as f:
        json.dump(sequence_lib_unique, f, indent=4)
    print(f"Unique sequences: {len(sequence_lib_unique)}")

    # 4.2 Write FASTA files with identical chain_A and chain_B entries
    fasta_output_path = Path("homodimer_fasta_lib")
    fasta_output_path.mkdir(parents=True, exist_ok=True)

    for seq_id, single_chain_seq in sequence_lib_unique.items():
        fasta_content = (
            f">protein|name=chain_A\n{single_chain_seq}\n"
            f">protein|name=chain_B\n{single_chain_seq}\n"
        )
        fasta_file_path = fasta_output_path / f"{seq_id}"
        with open(fasta_file_path, "w") as fasta_file:
            fasta_file.write(fasta_content)
        print(f"Written {fasta_file_path}")




def generate_redesign_heterodimer(mutation_combinations, redesign_sequence):
    
    wildtype_seq = redesign_sequence.copy()
    mutated_seq_lib= {}
    for i, mutation_set in enumerate(mutation_combinations, start=1):
        print(f"Processing mutation set {i}: {mutation_set}")

        # replace the residues in the wildtype sequence with the mutations
        for pos,new_aa in mutation_set.items():
            # update the wildtype sequence with mutations
            wildtype_seq[int(pos)] = new_aa
            
            # add the mutated sequence to the library
        mutated_seq_lib[f'mutated_sequence_{i}'] = "".join(wildtype_seq.copy().values())
        
    return mutated_seq_lib


    
def do_heterodimer_pipeline(
    pdb_name: str,
    pocket_size: float,
    cutoff: float,
    probabilities_path: str = "cnn_timed_50percentmodel_cofactor-148-1.73.csv",
    label_map_path: str = "datasetmap.txt",    
) -> Dict[str, str]:

    '''
    for heterodimer design

    Returns a dict of mutated sequence name -> sequence (unique if dedupe=True).
    Also writes:
      - {pdb_name}_predict_prob.csv
      - {pdb_name}_redesign_pocket_sequence.json
      - mutated_heterodimer_sequences.json

    '''
    
    # step1: Process generated CSV file from timed design CSV ---
    
    output_csv_name = f'{pdb_name}_predict_prob.csv'

    df_prob = process_data(probabilities_path, label_map_path, output_csv_name,pdb_name)

    # --- Section 2: Process the redesign sequence and JSON file from the selected sequence ---
    get_pocket_residues(pocket_size)
    current_directory = os.getcwd()

    # Load redesign sequence information from a JSON file

    json_input_path = os.path.join(current_directory, 'pocket_residues.json')
    with open(json_input_path, 'r') as f:
        redesign_sequence = json.load(f)
    predict_seq = {}
    for key_full_name in redesign_sequence.keys():

        predict_seq[key_full_name] = {}  # Initialize a sub-dictionary for each key_full_name
        for item_idx in redesign_sequence[key_full_name]:
            # Ensure item_idx is an integer for lookup
            predict_seq[key_full_name][item_idx] = get_amino_acid(key_full_name, int(item_idx), df_prob, cutoff)
    
    # Save the predicted sequence to a JSON file
    output_json_name = f'{pdb_name}_redesign_pocket_sequence.json'
    json_output_path = os.path.join(current_directory, output_json_name)
    
    with open(json_output_path, 'w') as f:
        json.dump(predict_seq, f, indent=4)
    
    print(f"Predicted sequence data saved to '{output_json_name}'.")

        
    # --- Section 3: Generate the mutatation combinations based on the predicted JSON file ---
    # Generate the mutation combinations
    process_pdb = pdb_name + '.pdb'
    redesign_sequence = generate_fasta(process_pdb)
    mutation_combinations = generate_combination(pdb_name, json_output_path)

    # --- Section 4: Generate the mutated heterodimer sequences based on the mutation combinations and
    # mutate back to the wildtype sequence
    mutated_all = generate_redesign_heterodimer(mutation_combinations, redesign_sequence)

    # check if there has repeat sequences
    existing_sequences = set()
    sequence_lib_unique = {k: v for k, v in mutated_all.items() if v not in existing_sequences and not existing_sequences.add(v)}


    # 4.1 convert the mutated_all to json file
    with open('mutated_heterodimer_sequences.json', 'w') as f:
        json.dump(sequence_lib_unique, f, indent=4)
    print(sequence_lib_unique)




    # 4.2 split sequence into chain A and chain B and generate a new fasta file
    fasta_output_path = 'heterodimer_fasta_lib'
    my_structure = ampal.load_pdb(process_pdb)
    chain_A_length = len(my_structure[0].sequence)
    chain_B_length = len(my_structure[1].sequence)

    print(f"Chain A length: {chain_A_length}, Chain B length: {chain_B_length}")

    for seq_id, full_sequence in sequence_lib_unique.items():
        chain_A_seq = full_sequence[:chain_A_length]
        chain_B_seq = full_sequence[chain_A_length:chain_A_length + chain_B_length]
        
        fasta_content = f">protein|name=chain_A\n{chain_A_seq}\n>protein|name=chain_B\n{chain_B_seq}\n"
        
        fasta_file_path = Path(fasta_output_path) / f"{seq_id}"
        # generate fasta path if not exist
        fasta_file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(fasta_file_path, 'w') as fasta_file:
            fasta_file.write(fasta_content)
        print(f"Written {fasta_file_path}")







def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Sequence redesign pipeline (heterodimer or homodimer)."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True, help="Choose pipeline")

    # Heterodimer subcommand
    p_hetero = subparsers.add_parser("heterodimer", help="Run heterodimer redesign")
    p_hetero.add_argument("--pdb_name", required=True, help="PDB name (without .pdb)")
    p_hetero.add_argument("--pocket_size", type=float, required=True, help="Pocket size (Å)")
    p_hetero.add_argument("--cutoff", type=float, required=True, help="Probability cutoff")
    p_hetero.add_argument("--prob_csv", default="cnn_timed_50percentmodel_cofactor-148-1.73.csv",
                          help="Predicted probabilities CSV from TIMED")
    p_hetero.add_argument("--label_map", default="datasetmap.txt", help="Label map file")
    
    
    # Homodimer subcommand
    p_homo = subparsers.add_parser("homodimer", help="Run homodimer redesign")
    p_homo.add_argument("--protein_name", required=True, help="Protein name (without .pdb)")
    p_homo.add_argument("--pocket_size", type=float, required=True, help="Pocket size (Å)")
    p_homo.add_argument("--cutoff", type=float, required=True, help="Probability cutoff")
    p_homo.add_argument("--prob_csv", default="cnn_timed_50percentmodel_cofactor-148-1.73.csv",
                        help="Predicted probabilities CSV from TIMED")
    p_homo.add_argument("--label_map", default="datasetmap.txt", help="Label map file")
    

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.mode == "heterodimer":
        result = do_heterodimer_pipeline(
            pdb_name=args.pdb_name,
            pocket_size=args.pocket_size,
            cutoff=args.cutoff,
            probabilities_path=args.prob_csv,
            label_map_path=args.label_map,
        )
        print("Processing heterodimer pipeline...")
    elif args.mode == "homodimer":
        result = do_homodimer_pipeline(
            protein_name=args.protein_name,
            pocket_size=args.pocket_size,
            cutoff=args.cutoff,
            probabilities_path=args.prob_csv,
            label_map_path=args.label_map,
        )
        print("Processing homodimer pipeline...")
    else:
        parser.error("Unknown mode.")



if __name__ == "__main__":
    main()