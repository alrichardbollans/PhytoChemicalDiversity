import os

import pandas as pd
from compile_trait_data import FAMILIES_OF_INTEREST
from wcvpy.wcvp_name_matching import get_accepted_info_from_names_in_column

from collect_and_compile_data.collect_compound_data import WCVP_VERSION


def prepare_duke_data():
    duke_file = os.path.join('inputs', 'DUKE_ETHNOBOT.csv')
    duke_df = pd.read_csv(duke_file)
    print(duke_df['ACTIVITY'].unique().tolist())
    duke_df = duke_df.dropna(subset=['ACTIVITY'])
    duke_df['full_name'] = duke_df['TAXON'] + " " + duke_df['TAXAUTHOR'].fillna('')
    accepted_duke = get_accepted_info_from_names_in_column(duke_df, 'full_name', wcvp_version=WCVP_VERSION)
    accepted_duke = accepted_duke[accepted_duke['accepted_family'].isin(FAMILIES_OF_INTEREST)]
    accepted_duke.to_csv('outputs/cleaned_DUKE_accepted.csv', index=False)


def count_number_of_medicinal_activities_in_region():
    accepted_duke = pd.read_csv('outputs/cleaned_DUKE_accepted.csv')
    region_data = pd.read_csv('../get_diversity_metrics/outputs/group_data/native_regions.csv', index_col=0)
    species_region_data = pd.read_csv('../get_diversity_metrics/outputs/group_info/native_regions.csv', index_col=0)
    out_data = []
    for region in region_data['Assigned_group'].unique():
        species_in_region = species_region_data[species_region_data['Assigned_group'] == region]
        medicinal_activities_in_region = accepted_duke[accepted_duke['accepted_species'].isin(species_in_region['accepted_species'].values)][
            'ACTIVITY'].unique().tolist()
        number_medicinal_activities_in_region = len(medicinal_activities_in_region)
        out_data.append([region, number_medicinal_activities_in_region])

    out_df = pd.DataFrame(out_data, columns=['Assigned_group', 'number_medicinal_activities_in_region'])
    out_df.to_csv('outputs/medicinal_activities_per_region.csv', index=False)


if __name__ == '__main__':
    prepare_duke_data()
    count_number_of_medicinal_activities_in_region()
