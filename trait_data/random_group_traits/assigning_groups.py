import os
from operator import index
from random import choice

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import hist
from phytochempy.chemical_diversity_metrics import get_pathway_based_diversity_measures, calculate_FAD_measures

from trait_data.collect_compound_data import resolve_compound_data_to_group, NP_PATHWAYS


def main():
    species_data = pd.read_csv(os.path.join('..', 'compile_trait_data', 'outputs', 'species_trait_data.csv'), index_col=0)[
        ['accepted_species', 'Genus', 'Animal Richness', 'PC0', 'PC1']]

    # For the moment, just mimic number of genera.
    number_of_groups = len(species_data['Genus'].unique().tolist())

    groups = range(0, number_of_groups)
    assigned_groups = [str(choice(groups)) for i in species_data.index]
    # hist(assigned_groups)
    # plt.show()
    species_data['Assigned_group'] = assigned_groups
    species_data = species_data.drop(columns=['Genus'])

    # resolve traits to group
    assert len(species_data[species_data.duplicated(subset=['accepted_species'])].index) == 0
    species_data['number_of_species_in_group'] = species_data[['Assigned_group', 'accepted_species']].groupby('Assigned_group').transform('count')

    mean_values = species_data[['Assigned_group', 'number_of_species_in_group', 'Animal Richness', 'PC0', 'PC1']].groupby(
        'Assigned_group').mean()
    mean_values = mean_values.reset_index()

    def check_means(x):
        if x != int(x):
            raise ValueError
        else:
            pass

    mean_values['number_of_species_in_group'].apply(check_means)
    print(mean_values)

    # After groups have been assigned, find compound data
    compound_data = pd.read_csv(os.path.join('..', 'collect_compound_data', 'outputs', 'all_species_compound_data.csv'), index_col=0)
    working_data = pd.merge(compound_data, species_data, how='left', on='accepted_species', validate='many_to_one')

    assert len(working_data) == len(compound_data)
    # then need to calculate metrics.
    group_compound_data, group_pathway_data = resolve_compound_data_to_group(working_data, 'Assigned_group')

    abundance_diversity = get_pathway_based_diversity_measures(group_pathway_data, NP_PATHWAYS, taxon_name_col='Assigned_group')
    FAD_measures = calculate_FAD_measures(group_compound_data, taxon_grouping='Assigned_group')

    compiled_data = pd.merge(mean_values, abundance_diversity, how='left', on='Assigned_group', validate='one_to_one')
    compiled_data = pd.merge(compiled_data, FAD_measures, how='left', on='Assigned_group', validate='one_to_one')

    return compiled_data


if __name__ == '__main__':
    main()
