import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

from diversity_metrics import genus_abundance_diversity_data_csv, genus_distance_diversity_data_csv, species_richness_csv


def main():
    ## Note this won't yet look the same as metric correlations as phy diversity isn't yet calculated for all genera
    distance_diversity_df = pd.read_csv(genus_distance_diversity_data_csv, index_col=0)
    abundance_diversity_df = pd.read_csv(genus_abundance_diversity_data_csv, index_col=0)
    richness_df = pd.read_csv(species_richness_csv, index_col=0)
    diversity_df = pd.merge(abundance_diversity_df, distance_diversity_df, on='Genus')
    diversity_df = pd.merge(diversity_df, richness_df, on='Genus')

    phylogenetic_measures = pd.read_csv(os.path.join('outputs', 'phylogenetic_diversities.csv'))[[
        'Genus', 'phylogenetic_diversity', 'genus_age', 'number_of_species_in_data_and_tree']]
    rename_dict = {'phylogenetic_diversity': 'PhyDiv', 'genus_age': 'GenusAge', 'number_of_species_in_data_and_tree': 'NumSp'}
    phylogenetic_measures = phylogenetic_measures.rename(columns=rename_dict)
    diversity_df = pd.merge(diversity_df, phylogenetic_measures, on='Genus')

    indices = ['FAD', 'MFAD', 'APWD', 'H', 'Hbc', 'J', 'G', 'PhyDiv', 'GenusAge', 'NumSp']
    ### Abundances
    corr_df = diversity_df[indices + ['N']].corr(method='kendall')
    corr_df.to_csv(os.path.join('outputs', 'correlations.csv'))
    # plot the heatmap

    corr_df = corr_df.drop(columns=['N'])
    corr_df = corr_df.drop('FAD')
    # Getting the Upper Triangle of the co-relation matrix
    mask = np.zeros_like(corr_df)
    mask[np.triu_indices_from(mask)] = True

    # Want diagonal elements as well
    mask[np.diag_indices_from(mask)] = False
    seaborn.heatmap(corr_df, cmap='coolwarm', annot=True, mask=mask, cbar=False)
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', 'diversity_heatmap.jpg'), dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
