import os

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from collect_compound_data import NP_PATHWAYS, genus_distinct_pathway_data_csv

_output_path = resource_filename(__name__, 'outputs')
genus_abundance_diversity_data_csv = os.path.join(_output_path, 'genus_level_pathway_diversity_information.csv')


def get_pathway_based_diversity_measures_for_genera(measure_df: pd.DataFrame, pathways: list) -> pd.DataFrame:
    ## Read data for all pathways into

    ### Begin with Shannon index
    measure_df['shannon'] = 0
    for pathway in pathways:
        measure_df[f'ln_mean_identified_as_{pathway}'] = np.log(measure_df[f'mean_identified_as_{pathway}']).replace(-np.inf, 0)
        measure_df['shannon'] = measure_df['shannon'] + measure_df[f'mean_identified_as_{pathway}'] * measure_df[
            f'ln_mean_identified_as_{pathway}']

    measure_df['shannon'] = -measure_df['shannon']
    ## Bias corrected shannon
    # From chao_nonparametric_2003, following beck_comparing_2010.
    # Note that there are updated metrics for calculating coverage e.g. chao_coveragebased_2012
    measure_df['number_singletons'] = 0
    for pathway in pathways:
        measure_df.loc[measure_df[f'identified_{pathway}_count'] == 1, 'number_singletons'] += 1
    measure_df['sample_coverage'] = 1 - (measure_df['number_singletons'] / measure_df['identified_compounds_count'])

    measure_df['bc_shannon'] = 0
    for pathway in pathways:
        ###New log
        measure_df[f'ln_Cmean_identified_as_{pathway}'] = np.log(
            measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage'])
        addition = ((measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage'] * measure_df[f'ln_Cmean_identified_as_{pathway}']) /
                    (1 - (1 - measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage']) ** measure_df[
                        'identified_compounds_count'])).fillna(0)
        measure_df['bc_shannon'] = measure_df['bc_shannon'] + addition
    measure_df['bc_shannon'] = -measure_df['bc_shannon']

    # Simpson index also refered to as gini-simpson
    # Used in e.g. corre_evaluation_2023
    measure_df['simpson'] = 0
    for pathway in pathways:
        measure_df['simpson'] = measure_df['simpson'] + (
                measure_df[f'mean_identified_as_{pathway}'] * measure_df[f'mean_identified_as_{pathway}'])
    measure_df['simpson'] = 1 - measure_df['simpson']

    measure_df['number_of_apparent_categories'] = 0
    for pathway in pathways:
        measure_df[f'binary_identified_as_{pathway}'] = 0
        measure_df.loc[measure_df[f'identified_{pathway}_count'] > 0, f'binary_identified_as_{pathway}'] = 1
        measure_df['number_of_apparent_categories'] = measure_df['number_of_apparent_categories'] + measure_df[
            f'binary_identified_as_{pathway}']

    ## Pielou index normalises by the 'richness' for the given genus, i.e. the number of different pathways present
    ## as discussed in corre_evaluation_2023.
    measure_df['pielou'] = measure_df['shannon'] / (np.log(measure_df['number_of_apparent_categories']))
    measure_df['pielou'] = measure_df['pielou'].fillna(measure_df['shannon'])


    measure_df['norm_bc_shannon'] = measure_df['bc_shannon']/measure_df['bc_shannon'].max()
    measure_df['norm_shannon'] = measure_df['shannon']/measure_df['shannon'].max()
    measure_df = measure_df[['Genus', 'shannon', 'bc_shannon', 'simpson', 'pielou', 'norm_bc_shannon', 'norm_shannon']]

    return measure_df


def main():
    g_df = pd.read_csv(genus_distinct_pathway_data_csv, index_col=0)
    out_df = get_pathway_based_diversity_measures_for_genera(g_df, NP_PATHWAYS)
    out_df.to_csv(genus_abundance_diversity_data_csv)


if __name__ == '__main__':
    main()
