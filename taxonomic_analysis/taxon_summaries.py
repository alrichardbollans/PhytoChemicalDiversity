import os

import pandas as pd
from pkg_resources import resource_filename

from trait_data.collect_compound_data import FAMILIES_OF_INTEREST, species_in_study_csv

_output_path = resource_filename(__name__, 'outputs')


def main():
    from wcvpy.wcvp_download import get_all_taxa, wcvp_columns

    all_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True, version='12')

    genus_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Genus']
    genus_df.describe(include='all').to_csv(os.path.join(_output_path, 'genus_summary.csv'))

    species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']
    species_df.describe(include='all').to_csv(os.path.join(_output_path, 'species_summary.csv'))

    species_in_study = pd.read_csv(species_in_study_csv)

    issues = species_in_study[~species_in_study['accepted_species'].isin(species_df['accepted_species'].values)]
    assert len(issues['accepted_species'].tolist()) == 0

    species_df = species_df[species_df['accepted_species'].isin(
        species_in_study['accepted_species'].values)]
    species_df.describe(include='all').to_csv(os.path.join(_output_path, 'species_in_study_summary.csv'))


if __name__ == '__main__':
    main()
