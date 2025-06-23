import os

import pandas as pd


def get_working_data():
    phy_div_data = pd.read_csv(
        os.path.join('..', 'collect_and_compile_data', 'get_phylogenetic_diversities', 'outputs', 'group_data', 'native_regions_transformed.csv'),
        index_col=0)
    trait_data = pd.read_csv(
        os.path.join('..', 'collect_and_compile_data', 'get_diversity_metrics', 'outputs', 'group_data', 'native_regions_transformed.csv'),
        index_col=0)

    trait_data = trait_data.rename(columns={
        'Assigned_group': 'Group'
    })

    working_data = pd.merge(phy_div_data, trait_data, on='Group', how='inner')
    if 'number_of_species_in_data_and_tree' in working_data.columns:
        issues = working_data[working_data['number_of_species_in_group'] != working_data['number_of_species_in_data_and_tree']]

        if len(issues) > 0:
            print(issues)
            raise ValueError
    working_data = working_data.rename(columns={'phylogenetic_diversity': 'PD', 'number_of_species_in_group': 'SR'})
    return working_data
