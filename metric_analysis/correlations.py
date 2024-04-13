import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

from diversity_metrics import genus_abundance_diversity_data_csv, genus_distance_diversity_data_csv, species_richness_csv


def both():
    distance_diversity_df = pd.read_csv(genus_distance_diversity_data_csv, index_col=0)
    abundance_diversity_df = pd.read_csv(genus_abundance_diversity_data_csv, index_col=0)
    richness_df = pd.read_csv(species_richness_csv, index_col=0)
    diversity_df = pd.merge(abundance_diversity_df, distance_diversity_df, on='Genus')
    diversity_df = pd.merge(diversity_df, richness_df, on='Genus')


    indices = ['FAD', 'MFAD', 'APWD'] + ['shannon','bc_shannon', 'pielou', 'simpson']
    unbound_indices = ['MFAD', 'bc_shannon', 'shannon']
    bound_indices = ['APWD', 'pielou', 'simpson']
    ### Abundances
    corr_df = diversity_df[indices + ['N']].corr()
    corr_df.to_csv(os.path.join('outputs', 'correlations.csv'))
    # plot the heatmap

    # Getting the Upper Triangle of the co-relation matrix
    matrix = np.triu(corr_df)
    seaborn.heatmap(corr_df, cmap='coolwarm', annot=True, mask=matrix)
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', 'diversity_heatmap.jpg'), dpi=300)
    plt.close()

    for i in unbound_indices:
        seaborn.regplot(data=diversity_df, y=i, x='N', label=i)

    plt.legend()
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join('outputs', 'unbound_metrics.jpg'), dpi=300)
    plt.close()
    for i in bound_indices:
        seaborn.regplot(data=diversity_df, y=i, x='N', label=i)
    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join('outputs', 'bound_metrics.jpg'), dpi=300)
    plt.close()

    for i in ['norm_FAD', 'norm_MFAD', 'APWD'] + ['norm_bc_shannon', 'pielou', 'norm_shannon', 'simpson']:
        seaborn.regplot(data=diversity_df, y=i, x='N', label=i)

    plt.legend()
    plt.ylim([-0.1, 1.1])
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join('outputs', 'all_metrics_with_minmaxing.jpg'), dpi=300)
    plt.close()


def main():
    # abundance_measures()
    # distance_measures()
    both()


if __name__ == '__main__':
    main()
