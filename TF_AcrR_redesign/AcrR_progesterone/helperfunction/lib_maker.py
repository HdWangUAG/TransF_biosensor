# for selected candidates in df_top_3_both, get the corresponding seq from 3cyl_TES_unique_sequences.csv
import pandas as pd
import os
from dnachisel import *
import numpy as np

def convert_aa_to_dna(seq):

    def weighted_random_codon(codons_with_weights):
        codons, weights = zip(*codons_with_weights)
        total = sum(weights)
        normalized_weights = [weight / total for weight in weights]
        return np.random.choice(codons, p=normalized_weights)

    
    '''
    Handing W 2025-04-03
    A function to convert amino acid sequence to DNA sequence.
    Args:
        seq: str, the amino acid sequence.
    Returns:
        dna_seq: str, the DNA sequence.
    '''
    codon_usage = {
        'A': [('GCT', 0.18), ('GCC', 0.26), ('GCA', 0.23), ('GCG', 0.33)],
        'C': [('TGT', 0.46), ('TGC', 0.54)],
        'D': [('GAT', 0.63), ('GAC', 0.37)],
        'E': [('GAA', 0.68), ('GAG', 0.32)],
        'F': [('TTT', 0.58), ('TTC', 0.42)],
        'G': [('GGT', 0.35), ('GGC', 0.37), ('GGA', 0.13), ('GGG', 0.15)],
        'H': [('CAT', 0.57), ('CAC', 0.43)],
        'I': [('ATT', 0.49), ('ATC', 0.39), ('ATA', 0.12)],
        'K': [('AAA', 0.78), ('AAG', 0.22)],
        'L': [('TTA', 0.13), ('TTG', 0.14), ('CTT', 0.11), ('CTC', 0.10), ('CTA', 0.04), ('CTG', 0.48)],
        'M': [('ATG', 1.0)],
        'N': [('AAT', 0.49), ('AAC', 0.51)],
        'P': [('CCT', 0.20), ('CCC', 0.12), ('CCA', 0.19), ('CCG', 0.49)],
        'Q': [('CAA', 0.70), ('CAG', 0.30)],
        'R': [('CGT', 0.36), ('CGC', 0.36), ('CGA', 0.06), ('CGG', 0.11), ('AGA', 0.06), ('AGG', 0.05)],
        'S': [('TCT', 0.17), ('TCC', 0.15), ('TCA', 0.14), ('TCG', 0.14), ('AGT', 0.16), ('AGC', 0.25)],
        'T': [('ACT', 0.25), ('ACC', 0.36), ('ACA', 0.16), ('ACG', 0.23)],
        'V': [('GTT', 0.28), ('GTC', 0.20), ('GTA', 0.17), ('GTG', 0.35)],
        'W': [('TGG', 1.0)],
        'Y': [('TAT', 0.59), ('TAC', 0.41)]}



    
    dna_seq = ''
    for amino_acid in seq:
        # convert amino acid to DNA
        try:

            dna_seq += weighted_random_codon(codon_usage[amino_acid])
        # check if the amino acid is in the codon usage dictionary
        except KeyError:
            print(f"Warning: {amino_acid} not found in codon usage dictionary. Using 'NNN' as placeholder.")


    return dna_seq

    
def remove_bsai(dna_sequence):
    '''
    Handing W 2025-04-03
    A function to remove BsaI site from the DNA sequence.
    Args:
        seq: str, the DNA sequence.
    Returns:
        seq: str, the DNA sequence without BsaI site.
    '''
    
    
    problem = DnaOptimizationProblem(
        sequence=dna_sequence,
        constraints=[
            AvoidPattern("BsaI_site"),
            
        ],
        )  # Note: always use a codon optimisation specification with EnforceTranslation
    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS

    print(problem.constraints_text_summary())
    print(problem.objectives_text_summary())

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)

    final_sequence = problem.sequence  # string
    return final_sequence

