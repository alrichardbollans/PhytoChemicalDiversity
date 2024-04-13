import os

from pkg_resources import resource_filename

from collect_compound_data import FAMILIES_OF_INTEREST

_output_path = resource_filename(__name__, 'outputs')



def main():
    from wcvp_download import get_all_taxa, wcvp_columns

    all_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True)

    genus_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Genus']
    genus_df.describe(include='all').to_csv(os.path.join(_output_path, 'genus_summary.csv'))

    species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']
    species_df.describe(include='all').to_csv(os.path.join(_output_path, 'species_summary.csv'))

if __name__ == '__main__':
    main()
