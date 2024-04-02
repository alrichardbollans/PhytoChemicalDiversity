import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

from library_info_and_data_import import processed_compounds_output_path


def abundance_measures():
    abundance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_pathway_diversity_information.csv'), index_col=0)
    abundance_indices = ['bc_shannon', 'pielou', 'shannon', 'simpson']

    ### Abundances

    # plot the heatmap
    corr_df = abundance_diversity_df[abundance_indices + ['identified_compounds_count']].corr()
    # Getting the Upper Triangle of the co-relation matrix
    matrix = np.triu(corr_df)
    seaborn.heatmap(corr_df, cmap='coolwarm', annot=True, mask=matrix)
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', 'abundance_heatmap.jpg'), dpi=300)
    plt.close()
    seaborn.regplot(data=abundance_diversity_df, y='shannon', x='identified_compounds_count', label='shannon')
    seaborn.regplot(data=abundance_diversity_df, y='bc_shannon', x='identified_compounds_count', label='bc_shannon')

    plt.legend()
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join('outputs', 'unbound_abundance_metrics.jpg'), dpi=300)

    plt.close()
    seaborn.regplot(data=abundance_diversity_df, y='pielou', x='identified_compounds_count', label='pielou')
    seaborn.regplot(data=abundance_diversity_df, y='simpson', x='identified_compounds_count', label='simpson')
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.savefig(os.path.join('outputs', 'bound_abundance_metrics.jpg'), dpi=300)
    plt.close()


def distance_measures():
    distance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_distance_diversity_information.csv'), index_col=0)
    indices = ['FAD', 'MFAD', 'APWD']

    # plot the heatmap
    corr_df = distance_diversity_df[indices + ['N']].corr()
    # Getting the Upper Triangle of the co-relation matrix
    matrix = np.triu(corr_df)
    seaborn.heatmap(corr_df, cmap='coolwarm', annot=True, mask=matrix)
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', 'distance_heatmap.jpg'), dpi=300)

    plt.close()

    seaborn.regplot(data=distance_diversity_df, y='MFAD', x='N', label='MFAD')
    seaborn.regplot(data=distance_diversity_df, y='FAD', x='N', label='FAD')
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.legend()
    plt.savefig(os.path.join('outputs', 'unbound_distance_metrics.jpg'), dpi=300)
    plt.close()

    seaborn.regplot(data=distance_diversity_df, y='APWD', x='N', label='APWD')
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.legend()
    plt.savefig(os.path.join('outputs', 'bound_distance_metrics.jpg'), dpi=300)
    plt.close()


def both():
    distance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_distance_diversity_information.csv'), index_col=0)
    abundance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_pathway_diversity_information.csv'), index_col=0)

    diversity_df = pd.merge(abundance_diversity_df, distance_diversity_df, on='Genus')
    indices = ['FAD', 'MFAD', 'APWD'] + ['bc_shannon', 'pielou', 'shannon', 'simpson']
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
        seaborn.regplot(data=diversity_df, y=i, x='identified_compounds_count', label=i)

    plt.legend()
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join('outputs', 'unbound_metrics.jpg'), dpi=300)
    plt.close()
    for i in bound_indices:
        seaborn.regplot(data=diversity_df, y=i, x='identified_compounds_count', label=i)
    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.xlabel('Number of Identified Compounds')
    plt.ylabel('Index')
    plt.savefig(os.path.join('outputs', 'bound_metrics.jpg'), dpi=300)
    plt.close()

    pd.testing.assert_series_equal(diversity_df['identified_compounds_count'], diversity_df['N'], check_names=False)


def main():
    abundance_measures()
    distance_measures()
    both()


if __name__ == '__main__':
    main()
