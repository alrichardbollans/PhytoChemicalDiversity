import os

import pandas as pd
from pkg_resources import resource_filename
from wcvpy.wcvp_download import plot_native_number_accepted_taxa_in_regions

from collect_and_compile_data.collect_compound_data import FAMILIES_OF_INTEREST, species_in_study_csv, WCVP_VERSION

_output_path = resource_filename(__name__, 'outputs')


def main():
    from wcvpy.wcvp_download import get_all_taxa, wcvp_columns

    all_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True, version=WCVP_VERSION)

    genus_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Genus']
    genus_df.describe(include='all').to_csv(os.path.join(_output_path, 'genus_summary.csv'))

    species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']
    species_df.describe(include='all').to_csv(os.path.join(_output_path, 'species_summary.csv'))

    species_in_study = pd.read_csv(species_in_study_csv)

    issues = species_in_study[~species_in_study['accepted_species'].isin(species_df['accepted_species'].values)]
    assert len(issues['accepted_species'].tolist()) == 0

    species_in_study = species_df[species_df['accepted_species'].isin(
        species_in_study['accepted_species'].values)]
    species_in_study.describe(include='all').to_csv(os.path.join(_output_path, 'species_in_study_summary.csv'))

    plot_native_number_accepted_taxa_in_regions(species_in_study, 'accepted_species', _output_path,
                                                'species_in_study_native_dist.jpg', wcvp_version=WCVP_VERSION)

    plot_native_number_accepted_taxa_in_regions(species_df, 'accepted_species', _output_path,
                                                'species_in_families_native_dist.jpg', wcvp_version=WCVP_VERSION)


if __name__ == '__main__':
    main()
