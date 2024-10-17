import os.path

from wcvpy.wcvp_download import get_all_taxa, get_distributions_for_accepted_taxa, wcvp_accepted_columns

from collect_and_compile_data.collect_compound_data import FAMILIES_OF_INTEREST, WCVP_VERSION
gentianales_species_distributions_csv = os.path.join('outputs', 'species_distributions.csv')

def main():
    acc_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True,
                            ranks=['Species'], version=WCVP_VERSION)
    wcvp_dists = get_distributions_for_accepted_taxa(acc_taxa,
                                                     wcvp_accepted_columns['name'],
                                                     include_extinct=True, wcvp_version=WCVP_VERSION).drop(
        columns=['accepted_parent_id', 'accepted_species_id']).sort_values(by=wcvp_accepted_columns['name']).reset_index(drop=True)
    wcvp_dists[['accepted_species','accepted_species_w_author','native_tdwg3_codes','intro_tdwg3_codes']].to_csv(gentianales_species_distributions_csv)


if __name__ == '__main__':
    main()
