import pandas as pd
from pkg_resources import resource_filename
from wcvp_download import wcvp_accepted_columns, plot_native_number_accepted_taxa_in_regions, get_all_taxa

from collect_compound_data import NP_PATHWAYS, FAMILIES_OF_INTEREST, processed_pathway_species_data_csv

_output_path = resource_filename(__name__, 'outputs')


def plot_underlying():
    family_data = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, ranks=['Species'], accepted=True)
    ## First plot the real underlying distribution
    plot_native_number_accepted_taxa_in_regions(family_data, wcvp_accepted_columns['species'], _output_path,
                                                'underlying_population.jpg', include_extinct=False)

    plot_native_number_accepted_taxa_in_regions(metabolite_data, wcvp_accepted_columns['species'], _output_path,
                                                'underlying_tested_population.jpg', include_extinct=False)


def plot_presences_of_compound(comp_class: str):
    # Then plot where class presences occur
    comp_examples = metabolite_data[metabolite_data[comp_class] == 1]
    plot_native_number_accepted_taxa_in_regions(comp_examples, wcvp_accepted_columns['species'], _output_path,
                                                comp_class + '.jpg', include_extinct=False)


if __name__ == '__main__':
    metabolite_data = pd.read_csv(processed_pathway_species_data_csv, index_col=0).dropna(subset=wcvp_accepted_columns['species'])
    plot_underlying()

    for comp_class in NP_PATHWAYS:
        plot_presences_of_compound(comp_class)
