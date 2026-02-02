__version__ = "1.0.0"
__author__ = "hdwang98"

import pandas as pd
import os
import glob
import matplotlib.pyplot as plt 
import seaborn as sns
import matplotlib.patches as mpatches

# calculate the ratio of holo to apo
def culling_structure(apo_score_path, holo_score_path, rosetta_energy_path, quality_threshold:  float = 0.8, confidence_threshold: float = 3.0):
    
    '''the function is designed to cull the protein structures based on their health scores by
    comparing the holo and apo scores, and then do ranking based on the sort criteria.
    '''

    sort_criteria = {
    
            # Cost effectiveness
            #'atp_cost_per_aa': -0.2,  
            #'sequence_length': -0.2,        
    
            # Codon variety  
            'dna_complexity_per_aa': 1,  
    
            # Detectability
            'sequence_molar_extinction_280':0.1,
    
            # esm metric
            #'pll_per_aa': 0.5,
    
            # boltz2 holo-protein health
            'complex_plddt_holo': 0.55,
            'mean_mut_plddt_holo': 1,
            'ptm_holo': -0.5,
            'complex_pde_holo': -0.25,

            # boltz2 apo-protein health
            'complex_plddt_apo': 0.55,
            'mean_mut_plddt_apo': 1,
            'ptm_apo': -0.5,      
            'complex_pde_apo': -0.25,
    
            # boltz2 holo/apo ratio
            'holo_apo_protein_ipde_ratio': -0.5, # lower is better
            'holo_apo_protein_iptm_ratio': 0.3,# higher is better
            'holo_apo_protein_iplddt_ratio': 0.3, # higher is better

            # holo_protein Rosetta energy
            'rosetta_energy_dG_dSASA_ratio':1.0,  # lower is better
    
    
            # # De Stress              
            # 'aggrescan3d_avg_value': -1,
            # 'hydrophobic_fitness': -0.5,
            # 'packing_density': 0.5,
            # 'rosetta_total_per_aa': -0.5,
    
            # NetSolP
            #'predicted_usability': 1,
        }

    # match the prot_name in holo and apo if exist
    def match_prot_name(apo_df, holo_df):

        matched_df = pd.merge(holo_df, apo_df, on='prot_name', suffixes=('_holo', '_apo'))
        
        return matched_df

    # load the csv files
    apo_df = pd.read_csv(apo_score_path)
    holo_df = pd.read_csv(holo_score_path)  
    rosetta_energy_df = pd.read_csv(rosetta_energy_path)
    # match the prot_name in holo and apo if exist
    matched_df = match_prot_name(apo_df, holo_df)



    merged_df_process = matched_df.copy()
    # calculate the health score based on the sort criteria for both holo and apo
    for key, value in sort_criteria.items():
        if key in matched_df.columns:
            if 'holo' in key:
                # extract the health score required items like plddt, pde, ptm
                merged_df_process[key] = matched_df[key] * value
            elif 'apo' in key:
                # extract the health score required items like plddt, pde, ptm
                merged_df_process[key] = matched_df[key] * value

    print("Processed DataFrame with Health Scores:")
    print(merged_df_process)

    merged_df_process['holo_apo_complex_ipde_ratio'] = merged_df_process['complex_ipde_holo'] / merged_df_process['complex_ipde_apo']
    merged_df_process['holo_apo_protein_iptm_ratio'] = merged_df_process['protein_iptm_holo'] / merged_df_process['protein_iptm_apo']
    merged_df_process['holo_apo_protein_iplddt_ratio'] = merged_df_process['complex_iplddt_holo'] / merged_df_process['complex_iplddt_apo']

    # 1. quality of holo and apo complex, for holo 
    merged_df_process['holo_complex_quality_score'] =(
        merged_df_process['complex_plddt_holo']  +
        merged_df_process['mean_mut_plddt_holo']  +
        merged_df_process['ptm_holo']  +
        merged_df_process['complex_pde_holo'] 

    )

    # 2. quality of holo and apo complex, for apo
    merged_df_process['apo_complex_quality_score'] =(
        merged_df_process['complex_plddt_apo'] +
        merged_df_process['mean_mut_plddt_apo']  +
        merged_df_process['ptm_apo']  +
        merged_df_process['complex_pde_apo']
    )

    # 3. calculate the Cofolding confidence score: -0.2*(holo_apo_complex_ipde_ratio) + 0.5*(holo_apo_protein_iptm_ratio) + 0.5*(holo_apo_protein_iplddt_ratio) 
    merged_df_process['Cofolding_confidence_score'] = (
        merged_df_process['holo_apo_complex_ipde_ratio'] +
        merged_df_process['holo_apo_protein_iptm_ratio'] +
        merged_df_process['holo_apo_protein_iplddt_ratio']

    )

    merged_df_process.to_csv('merged_apo_holo_boltz_scores.csv', index=False)


    df_high_quality_str = pd.DataFrame(columns=['prot_name', 'holo_complex_quality_score', 'apo_complex_quality_score', 'Cofolding_confidence_score'])
    merged_df_process = pd.read_csv('merged_apo_holo_boltz_scores.csv')

    # Filter the DataFrame for high-quality complexes, quality score larger than good score

    df_high_quality_str = merged_df_process[
        (merged_df_process['holo_complex_quality_score'] > quality_threshold) &
        (merged_df_process['apo_complex_quality_score'] > quality_threshold)]


    # further filtering the binding confidence score larger than confidence_threshold
    df_high_quality_str = df_high_quality_str[df_high_quality_str['Cofolding_confidence_score'] > confidence_threshold]
    #df_high_quality_str.to_csv('high_quality_complexes.csv', index=False)


    # again, we add the rosetta energy score to do the ranking
    ranking_df = pd.merge(df_high_quality_str, rosetta_energy_df, on='prot_name', how='left')
    ranking_df = ranking_df.sort_values(by=['Cofolding_confidence_score', 'dG_SASA_ratio'], ascending=[False, True])
    ranking_df.to_csv('ranked_rosetta_quality_complexes.csv', index=False)



    
    return merged_df_process

def plot_quality_distributions(merged_df_process):
    plt.figure(figsize=(10, 6))
    plt.hist(merged_df_process['holo_complex_quality_score'], bins=30, color='green', alpha=0.7, label='Holo Complex Quality')
    plt.hist(merged_df_process['apo_complex_quality_score'], bins=30, color='orange', alpha=0.7, label='Apo Complex Quality')
    plt.title('Distribution of Holo and Apo Complex Quality')
    plt.xlabel('Quality Score')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    #plt.savefig('complex_quality_distribution.png')




    plt.figure(figsize=(10, 6))
    plt.hist(merged_df_process['Cofolding_confidence_score'], bins=30, color='blue', alpha=0.7)
    plt.title('Distribution of Cofolding Confidence Scores')
    plt.xlabel('Cofolding Confidence Score')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)
    #plt.savefig('Cofolding_confidence_score_distribution.png')   

    plt.show()


# visualize the mutation distributation 
def count_mutation(df,wt_seq:str)->dict:
    
    '''
    the function is designed to count the mutations in a given DataFrame of protein sequences
    compared to a wildtype sequence.

    output: a dict
    {
        "prot_name": {
            "position": [(pos1, mutated_aa1), (pos2, mutated_aa2), ...],
            ...
        },
        ...
    }
    a figure shown in bar plot 

    revise at 13/08/25
    '''
    
    mutation_lib = {}

    for index, row in df.iterrows():
        prot_name = row['prot_name']
        sequence = str(row['sequence']) # Ensure it's a string

        min_len = min(len(wt_seq), len(sequence))
        if len(wt_seq) != len(sequence):
            print(f"Warning: Length mismatch for {prot_name}. WT: {len(wt_seq)}, Current: {len(sequence)}. Comparing up to shortest length.")

        mutations = [(i + 1, b) for i, (a, b) in enumerate(zip(wt_seq[:min_len], sequence[:min_len])) if a != b]

        # Save into the mutation_lib
        mutation_lib[prot_name] = mutations


    plot_data = []
    for prot_name, mutations in mutation_lib.items():
        for pos, mutated_aa in mutations:
            plot_data.append({'Position': pos, 'Mutated_AA': mutated_aa})

    # for same position, count the num of each type of mutation
    mutation_counts = {}
    for entry in plot_data:
        position = entry['Position']
        mutated_aa = entry['Mutated_AA']
        if position not in mutation_counts:
            mutation_counts[position] = {}
        if mutated_aa not in mutation_counts[position]:
            mutation_counts[position][mutated_aa] = 0
        mutation_counts[position][mutated_aa] += 1





    # plot 
    def plot_stack_bar(aggregated_mutation_data):
        plot_data = []
        for position, mutations_at_pos in aggregated_mutation_data.items():
            for mutated_aa, frequency in mutations_at_pos.items():
                plot_data.append({
                    'Position': position,
                    'Mutated_AA': mutated_aa,
                    'Frequency': frequency
                })

        df_plot = pd.DataFrame(plot_data)
        print(df_plot.head())
        df_plot['Position'] = df_plot['Position'].astype(int)
        df_plot = df_plot.sort_values(by='Position')

        plt.figure(figsize=(18, 7))
        # Use seaborn.barplot for evenly spaced categorical x-axis
        ax = sns.barplot(
            data=df_plot,
            x='Position',
            y='Frequency',
            hue='Mutated_AA',
            dodge=False,
        )
        

        # To stack the bars correctly with barplot, we need to manually plot each layer
        unique_positions = df_plot['Position'].unique()
        bottoms = {pos: 0 for pos in unique_positions}
        
        # Sort the dataframe to ensure consistent stacking order
        df_plot = df_plot.sort_values(by=['Position', 'Mutated_AA'])
        
        # Create a custom legend with a specific location to avoid overlap
        # Define the 20 standard amino acids in a fixed order
        standard_aas = [
            'A', 'R', 'N', 'D', 'C',
            'Q', 'E', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P',
            'S', 'T', 'W', 'Y', 'V'
        ]
        # Consistent color map
        palette = sns.color_palette("husl", n_colors=len(standard_aas))
        color_map = dict(zip(standard_aas, palette))
        ax = plt.gca() # Get the current axes
        
        # Iterate through each unique mutated amino acid to create a stacked layer
        for mutated_aa, group_df in df_plot.groupby('Mutated_AA'):
            sns.barplot(
                x='Position',
                y='Frequency',
                data=group_df,
                ax=ax,
                bottom=[bottoms[pos] for pos in group_df['Position']], # Use `bottom` parameter for stacking
                color=color_map[mutated_aa],
                label=mutated_aa
            )
            # Update the bottoms for the next layer
            for i, row in group_df.iterrows():
                bottoms[row['Position']] += row['Frequency']

        # --- Start of the fix for the blank space ---
        # Get the unique positions and create an index-based tick distribution
        positions = df_plot['Position'].unique()
        ax.set_xticks(range(len(positions))) # Set ticks at integer indices
        ax.set_xticklabels(positions, rotation=45, fontsize=14) # Set labels to the actual position numbers
        # --- End of the fix ---

        plt.title('Redesigned 2oys Stacked Mutation Frequency by Position', fontsize=24)
        plt.xlabel('Position', fontsize=20)
        plt.ylabel('Total Mutations', fontsize=20)
        plt.yticks(fontsize=16)

        # revise the legend
        # Manual legend
        patches = [mpatches.Patch(color=color_map[aa], label=aa) for aa in standard_aas]
        ax.legend(handles=patches, title='Mutated Amino Acids', loc='upper left', bbox_to_anchor=(1, 1))


        plt.tight_layout()
        plt.show()
    
    # Call the plotting function
    plot_stack_bar(mutation_counts)

    

    return mutation_counts







# Example usage
'''
apo_score_path = 'path/to/apo_scores.csv'
holo_score_path = 'path/to/holo_scores.csv'
# rescoring, use a new weights, and do culling
merged_df_process, df_high_quality_str = culling_structure(apo_score_path, holo_score_path)
plot_quality_distributions(merged_df_process)
wt_seq = "ACDEFGHIKLMNPQRSTVWY"  # Example wildtype sequence, replace with your actual wildtype sequence
count_mutation(df_high_quality_str, wt_seq) # replace with your actual wildtype sequence


'''
