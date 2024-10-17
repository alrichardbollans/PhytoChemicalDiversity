import itertools
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

from collect_and_compile_data.get_diversity_metrics import genus_distance_diversity_data_csv, genus_abundance_diversity_data_csv, species_richness_csv


def both(diversity_df, tag: str):
    out_dir = os.path.join('outputs', tag)
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    indices = ['FAD', 'MFAD', 'APWD'] + ['H', 'Hbc', 'J', 'G']
    unbound_indices = ['MFAD', 'Hbc', 'H']
    bound_indices = ['APWD', 'J', 'G']
    ### Abundances
    corr_df = diversity_df[indices + ['N']].corr()
    corr_df.to_csv(os.path.join(out_dir, 'correlations.csv'))
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
    plt.savefig(os.path.join(out_dir, 'diversity_heatmap.jpg'), dpi=300)
    plt.close()

    for i in unbound_indices:
        seaborn.regplot(data=diversity_df, y=i, x='N', label=i)

    plt.legend()
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join(out_dir, 'unbound_metrics.jpg'), dpi=300)
    plt.close()
    for i in bound_indices:
        seaborn.regplot(data=diversity_df, y=i, x='N', label=i)
    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join(out_dir, 'bound_metrics.jpg'), dpi=300)
    plt.close()

    for i in ['FAD_minmax', 'MFAD_minmax', 'APWD_minmax'] + ['Hbc_minmax', 'J_minmax', 'H_minmax', 'G_minmax']:
        seaborn.regplot(data=diversity_df, y=i, x='N', label=i)

    plt.legend()
    plt.ylim([-0.1, 1.1])
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join(out_dir, 'all_metrics_with_minmaxing.jpg'), dpi=300)
    plt.close()


def main():
    distance_diversity_df = pd.read_csv(genus_distance_diversity_data_csv, index_col=0)
    abundance_diversity_df = pd.read_csv(genus_abundance_diversity_data_csv, index_col=0)
    richness_df = pd.read_csv(species_richness_csv, index_col=0)
    diversity_df = pd.merge(abundance_diversity_df, distance_diversity_df, on='Genus')
    genus_diversity_df = pd.merge(diversity_df, richness_df, on='Genus')

    both(genus_diversity_df, 'genera')

    region_trait_data = pd.read_csv(os.path.join('..', 'collect_and_compile_data', 'other_group_traits', 'outputs', 'group_data', 'native_regions.csv'),
                                    index_col=0)
    both(region_trait_data, 'native_regions')



if __name__ == '__main__':
    main()
