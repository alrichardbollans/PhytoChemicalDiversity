import os.path
from typing import List

import pandas as pd
from pkg_resources import resource_filename

from trait_data.collect_compound_data import FAMILIES_OF_INTEREST, WCVP_VERSION, species_in_study_csv

ENVIRON_VARS = ['Bio1', 'Bio4', 'Bio10', 'Bio11', 'Bio12', 'Bio15', 'Bio16', 'Bio17', 'Brkl Elevation', 'Elevation', 'Slope', 'Soil Nitrogen',
                'Soil pH', 'Soil Depth', 'Soil OCS', 'Soil Water Cap', 'Latitude', 'Longitude']

_output_path = resource_filename(__name__, 'outputs')

_wcvp_species_for_families_csv = os.path.join(_output_path, 'wcvp_data_species_for_families.csv')
wcvp_species_data_for_study_species_csv = os.path.join(_output_path, 'wcvp_species_data_for_study_species.csv')


def merge_new_vars_from_data(in_df: pd.DataFrame, var_names: List[str], var_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge new variables from data into the given DataFrame.

    :param var_df:
    :param in_df: The input DataFrame to merge new variables into.
    :param var_names: A list of variable names to merge.
    :return: The updated DataFrame with the new variables merged.
    """

    for var in var_names:
        if var not in var_df.columns:
            print(var_df.columns)
            raise ValueError(var)

    # Check there's no species dups
    duplicates = var_df[var_df.duplicated([wcvp_accepted_columns['species']], keep=False)]
    if len(duplicates.index) > 0:
        # If there is, make sure only have species then recheck
        var_df = var_df[var_df[wcvp_accepted_columns['rank']] == 'Species']
        duplicates = var_df[var_df.duplicated([wcvp_accepted_columns['species']], keep=False)]
        if len(duplicates.index) > 0:
            raise ValueError

    # check fams
    var_df = var_df[var_df['accepted_family'].isin(FAMILIES_OF_INTEREST)]

    assert len(var_df['accepted_family'].dropna().unique().tolist()) == 5

    updated_df = pd.merge(in_df, var_df[[wcvp_accepted_columns['species']] + var_names], how='left',
                          left_on=wcvp_accepted_columns['species'],
                          right_on=wcvp_accepted_columns['species'])

    return updated_df


def main():
    from wcvpy.wcvp_download import get_all_taxa

    all_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True, version=WCVP_VERSION)

    updated_trait_df = all_taxa[all_taxa[wcvp_accepted_columns['rank']] == 'Species']

    updated_trait_df['Genus'] = updated_trait_df[wcvp_accepted_columns['parent_name']]

    ### Lifeforms
    lifeform_data = pd.read_csv(os.path.join('inputs', 'lifeforms.csv'))
    habit_cols = ['Herb', 'Liana', 'Succulent', 'Shrub', 'Subshrub', 'Tree']
    updated_trait_df = merge_new_vars_from_data(updated_trait_df, habit_cols,
                                                lifeform_data)

    ### Lifeforms
    animal_data = pd.read_csv(os.path.join('inputs', 'mean_animal_region_richness_for_plants.csv'))
    # for lf in LIFEFORM_COLS:
    updated_trait_df = merge_new_vars_from_data(updated_trait_df, ['Animal Richness'],
                                                animal_data)

    #### Species
    clim_data = pd.read_csv(os.path.join('inputs', 'compiled_climate_vars.csv'))

    updated_trait_df = merge_new_vars_from_data(updated_trait_df, ENVIRON_VARS,
                                                clim_data)

    # Aggregated environ vars
    data_from_known_regions = pd.read_csv(os.path.join('inputs', 'compiled_climate_vars_from_known_regions.csv'))[
        [wcvp_accepted_columns['species']] + ENVIRON_VARS].set_index(
        wcvp_accepted_columns['species'])
    updated_trait_df = updated_trait_df.set_index(wcvp_accepted_columns['species'])
    updated_trait_df.update(data_from_known_regions, overwrite=False)

    # TODO: Fit PCA on all data to get PC columns

    # restrict to study species
    species_in_study = pd.read_csv(species_in_study_csv)
    updated_trait_df = updated_trait_df[updated_trait_df['accepted_species'].isin(
        species_in_study['accepted_species'].values)]

    out_df = updated_trait_df.reset_index()

    out_df.to_csv(os.path.join('outputs', 'species_trait_data.csv'))

    mean_values = out_df[['Genus', 'Animal Richness'] + habit_cols + ENVIRON_VARS].groupby('Genus').mean()

    print(mean_values)

    mean_values.to_csv(os.path.join('outputs', 'genus_trait_data.csv'))


if __name__ == '__main__':
    from wcvpy.wcvp_download import get_all_taxa, wcvp_accepted_columns, wcvp_columns

    main()
