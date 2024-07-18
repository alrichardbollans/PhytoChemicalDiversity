import os

import pandas as pd
from phytochempy.chemical_diversity_metrics import calculate_FAD_measures
from pkg_resources import resource_filename

from collect_compound_data import all_genus_compound_csv

_output_path = resource_filename(__name__, 'outputs')
genus_distance_diversity_data_csv = os.path.join(_output_path, 'genus_level_distance_diversity_information.csv')
if not os.path.isdir(_output_path):
    os.mkdir(_output_path)


def main():
    my_df = pd.read_csv(all_genus_compound_csv, index_col=0)
    FAD_measures = calculate_FAD_measures(my_df, taxon_grouping='Genus')
    FAD_measures.to_csv(genus_distance_diversity_data_csv)


if __name__ == '__main__':
    main()
