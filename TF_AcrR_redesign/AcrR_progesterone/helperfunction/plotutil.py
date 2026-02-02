__version__ = "0.1.0"
# updating at 20250709

import seaborn as sns
import matplotlib.pyplot as plt 
import pandas as pd
from scipy.stats import pearsonr



def plot_pearsonr_scatter(df: pd.DataFrame, prot_name: str, x: str , y :str ):
    """
    Plot the correlation between dG_SASA_ratio and protein_iptm.
    
    Parameters:
        df (pd.DataFrame): The dataframe containing the data.
        prot_name (str): The name of the protein.
    """
    
    # check if the required columns are in the dataframe
    if x not in df.columns or y not in df.columns:
        raise ValueError(f"Dataframe must contain '{x}' and '{y}' columns.")


    plt.figure(figsize=(12, 8))
    # use regplot to plot the regression line
    sns.regplot(
        x= x, 
        y= y, 
        data=df, 
        scatter=False, 
        line_kws={'color': 'red', 'linestyle': '--'}
    )

    # use scatter to plot the points
    sc = plt.scatter(
        df[x], 
        df[y], 
        c=df[y], 
        cmap='viridis', 
        alpha=0.7
    )

    plt.xlabel(f'{x}', fontsize=18)
    plt.ylabel(f'{y}', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.grid(True)

    # add colorbar adn fontsize of legend
    cbar = plt.colorbar(sc, label=f'{y}')
    cbar.set_label(f'{y}', fontsize=16)
    # calculate the pearsonr
    r, p_value = pearsonr(df[x], df[y])

    plt.title(
        f'Correlation between {x} and {y} for {prot_name}\n'
        f'Pearson r = {r:.3f}, p-value = {p_value:.3e}', 
        fontsize=22
    )
    plt.tight_layout()
    plt.savefig(f'{prot_name}{x}_and_{y}.png')
    plt.show()

