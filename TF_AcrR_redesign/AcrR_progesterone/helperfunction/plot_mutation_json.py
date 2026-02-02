import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns # Seaborn is often used with Matplotlib for enhanced heatmaps

# for predicted mutations json files.
#  x axis is position, y is amino acids and color shows values

with open('/home/hdwang/protein_sensor_hd/1flm/eng_interface/1flm_HCY4_250/1flm_HCY4_250_predicted_mutations.json', 'r') as data:
    json_data = json.load(data)


# set key as x-axis, y is 20 amino acids type, value is ddG
all_amino_acids = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
]
# Get all unique positions from the input data
all_positions = list(json_data.keys())

# Create a dictionary to store ddG values, mapping (position, amino_acid) to ddG
ddg_map = {}
for position, mutations in json_data.items():
    for amino_acid, data in mutations.items():
        ddg_map[(position, amino_acid)] = data['ddG']

# Create a 2D numpy array for the heatmap data
# Rows will be amino acids, columns will be positions
heatmap_data = np.zeros((len(all_amino_acids), len(all_positions)))
heatmap_data[:] = np.nan # Initialize with NaN for blank cells

for i, aa in enumerate(all_amino_acids):
    for j, pos in enumerate(all_positions):
        if (pos, aa) in ddg_map:
            heatmap_data[i, j] = ddg_map[(pos, aa)]

# Plotting the heatmap using seaborn (which is built on matplotlib)
plt.figure(figsize=(len(all_positions) * 0.8 + 2, len(all_amino_acids) * 0.4 + 2)) # Adjust figure size dynamically
ax = sns.heatmap(
    heatmap_data,
    xticklabels=all_positions,
    yticklabels=all_amino_acids,
    cmap='viridis',
    cbar=True, # Ensure colorbar is drawn
    cbar_kws={'label': 'ddG Value'},
    linewidths=.5,
    linecolor='gray',
    mask=np.isnan(heatmap_data)
)

# Get the colorbar object
cbar = ax.collections[0].colorbar

# Set the fontsize of the colorbar label
cbar.set_label('ddG Value', fontsize=14) # You can adjust '14' to your desired size
ax.patch.set_edgecolor('black')  # Set the border color
ax.patch.set_linewidth(1.5)
# You might also want to increase the tick label size on the colorbar
cbar.ax.tick_params(labelsize=12) # Adjust '12' to your desired size for ticks

plt.title('Mutations scanning for 1flm_HCY4', fontsize=20)
plt.xlabel('Residue Position', fontsize=18)
plt.ylabel('Amino Acid Type', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=16)
# set the colorbar label font size

plt.tight_layout() # Adjust layout to prevent labels from overlapping
plt.show()
