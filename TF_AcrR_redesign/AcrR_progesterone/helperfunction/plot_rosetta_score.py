# Before filtering, we can plot the distribution of the scores
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
prot_name = "1flm"  # specify the protein name
# Load the data
df = pd.read_csv(f'{prot_name}_confidence_rosetta_results.csv')

# Plot the distribution of scores using boxplot

plt.figure(figsize=(12, 6))
sns.boxplot(data=df[[ 'confidence_score','ptm', 'iptm', 'ligand_iptm', 'protein_iptm',
                     'complex_plddt', 'complex_iplddt', 'complex_pde', 'complex_ipde']])

plt.title('Distribution of boltz-2.0 scores for ' + prot_name, fontsize=16)
yticks = [ 0, 0.5, 1, 1.5, 2, 2.5, 3,3.5,4,4.5,5,5.5, 6,6.5,7]
# hide outliers
plt.ylim(0,7)
plt.ylabel('Score Value',fontsize=14)
#plt.xlabel('Score Type')

# show the median value for each score type
medians = df[[ 'confidence_score','ptm', 'iptm', 'ligand_iptm', 'protein_iptm',
                     'complex_plddt', 'complex_iplddt', 'complex_pde']].median()
max_values = df[['confidence_score','ptm', 'iptm', 'ligand_iptm', 'protein_iptm',
                 'complex_plddt', 'complex_iplddt', 'complex_pde']].max()

for i, (median, max_v) in enumerate(zip(medians, max_values)):
    plt.text(
        i, 
        max_v + 0.1, 
        f'{median:.2f}', 
        ha='center', 
        color='black',
        fontsize=12
    )

# for complex_ipde, we  put median value on the center of the box
plt.text(
    8, 
    df['complex_ipde'].median() + 0.1,
    f'{df["complex_ipde"].median():.2f}',
    ha='center', 
    color='black',
    fontsize=12
)


plt.grid(axis='y')
plt.yticks(yticks, fontsize=12)
plt.xticks(rotation=45,fontsize=14)
plt.tight_layout()
plt.savefig(f'{prot_name}_score_distribution.png')
plt.show()
