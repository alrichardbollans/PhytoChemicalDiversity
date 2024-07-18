import os

import numpy as np
import pandas as pd
from phytochempy.chemical_diversity_metrics import get_pathway_based_diversity_measures
from pkg_resources import resource_filename

from collect_compound_data import NP_PATHWAYS, genus_distinct_pathway_data_csv

_output_path = resource_filename(__name__, 'outputs')
genus_abundance_diversity_data_csv = os.path.join(_output_path, 'genus_level_pathway_diversity_information.csv')


def main():
    g_df = pd.read_csv(genus_distinct_pathway_data_csv, index_col=0)
    out_df = get_pathway_based_diversity_measures(g_df, NP_PATHWAYS)
    out_df.to_csv(genus_abundance_diversity_data_csv)


if __name__ == '__main__':
    main()
