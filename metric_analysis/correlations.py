import itertools
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn

from library_info_and_data_import import processed_compounds_output_path


def your_function():
    diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, f'genus_level_pathway_diversity_information.csv'))
    indices = ['e_bias_corrected_shannon_index', 'bias_corrected_shannon_index', 'pielou_index', 'shannon_index', 'simpson_index',
               'normalised_simpson', 'normalised_shannon', 'identified_compounds_count']
    for pair in itertools.combinations(indices, 2):
        i = pair[0]
        j = pair[1]
        if i != j:
            correlation = diversity_df[[i, j]].corr().iloc[0, 1]
            print(f"Correlation between {i} and {j}: {correlation}")
    seaborn.regplot(data=diversity_df, y='shannon_index', x='identified_compounds_count', label='shannon')
    seaborn.regplot(data=diversity_df, y='bias_corrected_shannon_index', x='identified_compounds_count', label='bias_corrected_shannon')

    plt.legend()
    plt.show()
    plt.close()
    seaborn.regplot(data=diversity_df, y='pielou_index', x='identified_compounds_count', label='pielou_index')
    seaborn.regplot(data=diversity_df, y='normalised_shannon', x='identified_compounds_count', label='normalised_shannon')
    seaborn.regplot(data=diversity_df, y='simpson_index', x='identified_compounds_count', label='simpson_index')
    seaborn.regplot(data=diversity_df, y='normalised_simpson', x='identified_compounds_count', label='normalised_simpson')

    plt.ylim([-0.1, 1.1])
    plt.legend()
    plt.show()


def main():
    your_function()


if __name__ == '__main__':
    main()
