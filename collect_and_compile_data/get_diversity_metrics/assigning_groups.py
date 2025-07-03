import ast
import os
from random import choice, sample

import numpy as np
import pandas as pd

from collect_and_compile_data.collect_compound_data import species_in_study_csv
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import resolve_traits_to_group


def get_region_groups():
    working_data = species_data_with_dists.copy()

    working_data = working_data.explode('native_tdwg3_codes')
    print(working_data)
    working_data = working_data.rename(columns={'native_tdwg3_codes': 'Assigned_group'}).drop(columns=['Genus'])
    ## Just a check
    assert len(working_data['Assigned_group'].unique()) == number_of_native_regions

    resolve_traits_to_group(
        working_data,
        tag='native_regions')

def get_region_groups_only_medicinal_species():
    working_data = species_data_with_dists.copy()

    medicinal_species = pd.read_csv(os.path.join('..', 'compile_trait_data', 'outputs', 'cleaned_DUKE_accepted.csv'), index_col=0)
    working_data = working_data[working_data['accepted_species'].isin(medicinal_species['accepted_species'].values)]
    working_data = working_data.explode('native_tdwg3_codes')
    print(working_data)
    working_data = working_data.rename(columns={'native_tdwg3_codes': 'Assigned_group'}).drop(columns=['Genus'])
    ## Just a check
    assert len(working_data['Assigned_group'].unique()) < number_of_native_regions

    resolve_traits_to_group(
        working_data,
        tag='native_regions_medicinal_species')

def get_species_groups():
    working_data = species_data.copy()

    working_data = working_data.drop(columns=['Genus'])
    working_data['Assigned_group'] = working_data['accepted_species']
    resolve_traits_to_group(
        working_data,
        tag='accepted_species')

def write_random_group(number_of_groups, largest_group_size, tag: str):
    """
    Generates and assigns random groups to species data, then resolves traits for each group.

    This function creates a specified number of random groups of species, assigns these
    groups to species in the dataset, and processes the assignment to resolve traits
    based on the provided tag. The input dataset is modified and used for further
    processing.

    Arguments:
        number_of_groups: int
            The total number of groups to be generated.
        largest_group_size: int
            The upper limit for the size of each group.
        tag: str
            A tag identifier used to resolve traits for each group.

    Returns:
        None
    """

    possible_group_sizes = range(1, largest_group_size)

    groups = {}
    for i in range(0, number_of_groups):
        group_size = choice(possible_group_sizes)
        groups[str(i)] = sample(species_data['accepted_species'].unique().tolist(), group_size)

    working_data = species_data.copy()
    working_data['Assigned_group'] = working_data['accepted_species'].apply(lambda x: [c for c in groups if x in groups[c]])
    working_data = working_data.explode('Assigned_group')

    working_data = working_data.drop(columns=['Genus'])

    resolve_traits_to_group(
        working_data,
        tag=tag)


def main():
    # get_region_groups()
    # get_species_groups()
    # Mimic number of genera.
    # number_of_groups = len(species_data['Genus'].unique().tolist())
    # largest_group_size = species_data[['Genus', 'accepted_species']].groupby('Genus').transform('count').max().iloc[0]
    # write_random_group(number_of_groups, largest_group_size, tag='random_genera')

    ## Mimic number of regions
    count_df = species_data_with_dists.explode('native_tdwg3_codes')
    largest_region_count = count_df[['native_tdwg3_codes', 'accepted_species']].groupby('native_tdwg3_codes').transform('count').max().iloc[0]
    write_random_group(number_of_native_regions, largest_region_count, tag='random_regions')


if __name__ == '__main__':
    species_data = pd.read_csv(species_in_study_csv, index_col=0)[
        ['accepted_species', 'Genus']]

    distribution_data = pd.read_csv(os.path.join('..', 'compile_trait_data', 'outputs', 'species_distributions.csv'), index_col=0)
    distribution_data = distribution_data.dropna(subset=['native_tdwg3_codes'])[['accepted_species', 'native_tdwg3_codes']]
    distribution_data['native_tdwg3_codes'] = distribution_data['native_tdwg3_codes'].apply(lambda x: ast.literal_eval(x))

    species_data_with_dists = pd.merge(distribution_data, species_data, how='right', on='accepted_species', validate='one_to_one')
    ## Just a check
    all_native_regions = []
    for regions in species_data_with_dists['native_tdwg3_codes'].unique():
        for r in regions:
            all_native_regions.append(r)
    all_native_regions = list(set(all_native_regions))
    number_of_native_regions = len(all_native_regions)
    with np.errstate(divide='ignore', invalid='ignore'):
        main()
