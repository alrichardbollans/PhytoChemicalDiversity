import pandas as pd
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_accepted_columns, plot_native_number_accepted_taxa_in_regions, get_all_taxa

from trait_data.collect_compound_data import NP_PATHWAYS, FAMILIES_OF_INTEREST, all_genus_compound_csv

_output_path = resource_filename(__name__, 'outputs')


def plot_underlying():
    family_data = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, ranks=['Genus'], accepted=True, version='12')
    ## First plot the real underlying distribution
    plot_native_number_accepted_taxa_in_regions(family_data, wcvp_accepted_columns['name'], _output_path,
                                                'underlying_population.jpg', include_extinct=False, wcvp_version='12')

    plot_native_number_accepted_taxa_in_regions(metabolite_data, 'Genus', _output_path,
                                                'underlying_tested_population.jpg', include_extinct=False, wcvp_version='12')


def plot_presences_of_compound(comp_class: str):
    # Then plot where class presences occur
    comp_examples = metabolite_data[metabolite_data[comp_class] == 1]
    plot_native_number_accepted_taxa_in_regions(comp_examples, 'Genus', _output_path,
                                                comp_class + '.jpg', include_extinct=False, wcvp_version='12')


if __name__ == '__main__':
    # Note with current set up this uses data where genera with single compounds (N=1) have been removed.
    metabolite_data = pd.read_csv(all_genus_compound_csv, index_col=0)
    plot_underlying()

    for comp_class in NP_PATHWAYS:
        plot_presences_of_compound(comp_class)
