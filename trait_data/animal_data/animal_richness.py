import os.path
import pandas as pd
import numpy as np
from wcvpy.wcvp_download import get_all_taxa, get_distributions_for_accepted_taxa, wcvp_accepted_columns

from collect_compound_data import FAMILIES_OF_INTEREST

_distributions_csv = os.path.join('outputs', 'distributions.csv')
mean_animal_region_richness_for_plants_csv = os.path.join('outputs', 'mean_animal_region_richness_for_plants.csv')


def get_dists():
    version = '12'

    acc_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True,
                            ranks=['Genus'], version=version)
    wcvp_dists = get_distributions_for_accepted_taxa(acc_taxa,
                                                     wcvp_accepted_columns['name'],
                                                     include_extinct=True, wcvp_version=version).drop(
        columns=['accepted_parent_id', 'accepted_species_id']).sort_values(by=wcvp_accepted_columns['name']).reset_index(drop=True)
    wcvp_dists.to_csv(_distributions_csv)


def get_mean_richness_vals_for_genera():
    """

    Calculate the mean richness values for each genus based on animal region richness data.

    Returns:
        DataFrame: DataFrame containing the mean richness values for each genus.

    """
    import ast

    from wcvpy.wcvp_download import native_code_column

    # From ApmTraits v1.12.2
    region_richness_df = pd.read_csv(os.path.join('inputs', 'animal_region_richness.csv'), index_col=0).reset_index(drop=True)
    _outputted_dist_df = pd.read_csv(_distributions_csv, index_col=0).reset_index(drop=True)

    def reformat_dist_col(given_val):
        if given_val == given_val:
            out = list(ast.literal_eval(given_val))
        else:
            out = np.nan
        return out

    _outputted_dist_df[native_code_column] = _outputted_dist_df[native_code_column].apply(reformat_dist_col)

    multilabels = _outputted_dist_df[native_code_column].str.join('|').str.get_dummies()

    for region in multilabels.columns:
        region_value = \
            region_richness_df[region_richness_df['tdwg3_codes'] == region]['animal_richness'].iloc[0]
        multilabels[region] = multilabels[region] * region_value
    multilabels = multilabels.replace(0, np.NaN)
    all_regions = multilabels.columns
    _outputted_dist_df['animal richness'] = multilabels[all_regions].mean(axis=1)
    assert _outputted_dist_df['accepted_rank'].unique().tolist() == ['Genus']

    out_cols = ['accepted_family', 'accepted_name']
    _outputted_dist_df = _outputted_dist_df[out_cols + ['animal richness']]
    _outputted_dist_df = _outputted_dist_df.rename(columns={'animal richness': 'Animal_Richness'})
    _outputted_dist_df = _outputted_dist_df.rename(columns={'accepted_name':'Genus'})
    _outputted_dist_df.to_csv(mean_animal_region_richness_for_plants_csv)
    return _outputted_dist_df


if __name__ == '__main__':
    # get_dists()
    get_mean_richness_vals_for_genera()
