import os

import pandas as pd
from phytochempy.chemical_diversity_metrics import get_pathway_based_diversity_measures, calculate_FAD_measures
from pkg_resources import resource_filename

from collect_and_compile_data.collect_compound_data import NP_PATHWAYS, all_species_compound_csv, \
    resolve_compound_data_to_group

_output_path = resource_filename(__name__, 'outputs')
genus_diversity_data_csv = os.path.join(_output_path, 'genera', 'genera_diversities.csv')


def genera():
    genus_abundance_diversity_data_csv = os.path.join(_output_path, 'genera', 'genus_level_pathway_diversity_information.csv')
    genus_distance_diversity_data_csv = os.path.join(_output_path, 'genera', 'genus_level_distance_diversity_information.csv')
    my_df = pd.read_csv(all_species_compound_csv, index_col=0)

    group_compound_data, group_pathway_data = resolve_compound_data_to_group(my_df, 'Genus')

    abundance_diversity = get_pathway_based_diversity_measures(group_pathway_data, NP_PATHWAYS, taxon_grouping='Genus')
    abundance_diversity.to_csv(genus_abundance_diversity_data_csv)

    FAD_measures = calculate_FAD_measures(group_compound_data, taxon_grouping='Genus')
    FAD_measures.to_csv(genus_distance_diversity_data_csv)

    compiled_data = pd.merge(abundance_diversity, FAD_measures, how='left', on='Genus', validate='one_to_one')

    for g in list(set(group_compound_data['Genus'].values.tolist())) + list(set(FAD_measures['Genus'].values.tolist())):
        assert g in compiled_data['Genus'].values

    compiled_data.to_csv(genus_diversity_data_csv)
    return compiled_data


if __name__ == '__main__':
    genera()
