import os

from pkg_resources import resource_filename

from collect_compound_data import FAMILIES_OF_INTEREST

_inputs_path = resource_filename(__name__, 'inputs')


def main():
    from wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns

    all_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True)
    # species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']
    # species_df = species_df[[wcvp_columns['name'], wcvp_accepted_columns['name'], 'genus']]

    # species_list = species_df[wcvp_columns['name']].unique().tolist()

    # with open(os.path.join('inputs', 'species_list.txt'), 'w') as f:
    #     for line in species_list:
    #         f.write(f"{line}\n")

    genus_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Genus']

    genus_list = genus_df[wcvp_columns['name']].unique().tolist()
    with open(os.path.join('inputs', 'genus_list.txt'), 'w') as f:
        for line in genus_list:
            f.write(f"{line}\n")

    genus_df[[wcvp_columns['name'], wcvp_columns['family']]].to_csv(os.path.join('inputs', 'genus_family_list.csv'))

if __name__ == '__main__':
    main()
