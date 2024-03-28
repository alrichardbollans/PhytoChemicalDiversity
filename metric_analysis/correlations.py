import itertools
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

from library_info_and_data_import import processed_compounds_output_path


def abundance_measures():
    abundance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_pathway_diversity_information.csv'), index_col=0)
    abundance_indices = ['bias_corrected_shannon_index', 'pielou_index', 'shannon_index', 'simpson_index',
                         'normalised_simpson', 'normalised_shannon']

    ### Abundances
    N_col = ['identified_compounds_count']
    # for pair in itertools.combinations(abundance_indices + N_col, 2):
    #     i = pair[0]
    #     j = pair[1]
    #     if i != j:
    #         correlation = abundance_diversity_df[[i, j]].corr().iloc[0, 1]
    #         print(f"Correlation between {i} and {j}: {correlation}")

    # plot the heatmap
    corr_df = abundance_diversity_df[abundance_indices].corr()
    # Getting the Upper Triangle of the co-relation matrix
    matrix = np.triu(corr_df)
    seaborn.heatmap(corr_df, cmap='coolwarm', annot=True, mask=matrix)
    plt.tight_layout()
    plt.show()
    plt.close()
    seaborn.regplot(data=abundance_diversity_df, y='shannon_index', x='identified_compounds_count', label='shannon')
    seaborn.regplot(data=abundance_diversity_df, y='bias_corrected_shannon_index', x='identified_compounds_count', label='bias_corrected_shannon')

    plt.legend()
    plt.show()
    plt.close()
    seaborn.regplot(data=abundance_diversity_df, y='pielou_index', x='identified_compounds_count', label='pielou_index')
    seaborn.regplot(data=abundance_diversity_df, y='normalised_shannon', x='identified_compounds_count', label='normalised_shannon')
    seaborn.regplot(data=abundance_diversity_df, y='simpson_index', x='identified_compounds_count', label='simpson_index')
    seaborn.regplot(data=abundance_diversity_df, y='normalised_simpson', x='identified_compounds_count', label='normalised_simpson')

    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.show()


def distance_measures():
    distance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_distance_diversity_information.csv'), index_col=0)
    indices = ['FAD', 'MFAD', 'APWD']
    N_col = ['N']

    # plot the heatmap
    corr_df = distance_diversity_df[indices].corr()
    # Getting the Upper Triangle of the co-relation matrix
    matrix = np.triu(corr_df)
    seaborn.heatmap(corr_df, cmap='coolwarm', annot=True, mask=matrix)
    plt.tight_layout()
    plt.show()
    plt.close()

    seaborn.regplot(data=distance_diversity_df, y='MFAD', x='N', label='MFAD')
    seaborn.regplot(data=distance_diversity_df, y='APWD', x='N', label='APWD')

    plt.legend()
    plt.show()
    plt.close()

    seaborn.regplot(data=distance_diversity_df, y='APWD', x='N', label='APWD')

    plt.legend()
    plt.show()
    plt.close()


def both():
    distance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_distance_diversity_information.csv'), index_col=0)
    abundance_diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, 'genus_level_pathway_diversity_information.csv'), index_col=0)

    diversity_df = pd.merge(abundance_diversity_df, distance_diversity_df, on='Genus')
    indices= ['FAD', 'MFAD', 'APWD'] +['bias_corrected_shannon_index', 'pielou_index', 'shannon_index', 'simpson_index',
                         'normalised_simpson', 'normalised_shannon']
    unbound_indices = ['MFAD', 'bias_corrected_shannon_index','shannon_index']
    bound_indices = ['APWD', 'pielou_index', 'simpson_index',
                         'normalised_simpson', 'normalised_shannon']
    ### Abundances
    N_col = ['identified_compounds_count']
    corr_df = diversity_df[indices + N_col + ['N']].corr()
    corr_df.to_csv(os.path.join('outputs', 'correlations.csv'))
    # plot the heatmap
    corr_df = diversity_df[indices].corr()
    # Getting the Upper Triangle of the co-relation matrix
    matrix = np.triu(corr_df)
    seaborn.heatmap(corr_df, cmap='coolwarm', annot=True, mask=matrix)
    plt.tight_layout()
    plt.show()
    plt.close()

    for i in unbound_indices:
        seaborn.regplot(data=diversity_df, y=i, x='identified_compounds_count', label=i)

    plt.legend()
    plt.show()
    plt.close()
    for i in bound_indices:
        seaborn.regplot(data=diversity_df, y=i, x='identified_compounds_count', label=i)
    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.show()
    plt.close()


    pd.testing.assert_series_equal(diversity_df['identified_compounds_count'], diversity_df['N'], check_names=False)


def main():
    # abundance_measures()
    # distance_measures()
    both()

if __name__ == '__main__':
    main()
