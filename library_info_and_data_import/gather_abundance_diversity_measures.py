import os

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from library_info_and_data_import import get_pathway_categories

_inputs_path = resource_filename(__name__, 'inputs')

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')

processed_compounds_output_path = resource_filename(__name__, 'outputs')


def get_pathway_based_diversity_measures_for_genera(measure_df: pd.DataFrame, pathways: list):
    ## Read data for all pathways into

    measure_df = measure_df[measure_df['identified_compounds_count'] > 1]

    ### Begin with Shannon index
    measure_df['shannon_index'] = 0
    for pathway in pathways:
        measure_df[f'ln_mean_identified_as_{pathway}'] = np.log(measure_df[f'mean_identified_as_{pathway}']).replace(-np.inf, 0)
        measure_df['shannon_index'] = measure_df['shannon_index'] + measure_df[f'mean_identified_as_{pathway}'] * measure_df[
            f'ln_mean_identified_as_{pathway}']

    measure_df['shannon_index'] = -measure_df['shannon_index']
    ## Bias corrected shannon
    # From chao_nonparametric_2003, following beck_comparing_2010.
    # Note that there are updated metrics for calculating coverage e.g. chao_coveragebased_2012
    measure_df['number_singletons'] = 0
    for pathway in pathways:
        measure_df.loc[measure_df[f'identified_{pathway}_count'] == 1, 'number_singletons'] += 1
    measure_df['sample_coverage'] = 1 - (measure_df['number_singletons'] / measure_df['identified_compounds_count'])

    measure_df['bias_corrected_shannon_index'] = 0
    for pathway in pathways:
        ###New log
        measure_df[f'ln_Cmean_identified_as_{pathway}'] = np.log(
            measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage'])
        addition = ((measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage'] * measure_df[f'ln_Cmean_identified_as_{pathway}']) /
                    (1 - (1 - measure_df[f'mean_identified_as_{pathway}'] * measure_df['sample_coverage']) ** measure_df[
                        'identified_compounds_count'])).fillna(0)
        measure_df['bias_corrected_shannon_index'] = measure_df['bias_corrected_shannon_index'] + addition
    measure_df['bias_corrected_shannon_index'] = -measure_df['bias_corrected_shannon_index']

    # Simpson index also refered to as gini-simpson
    # Used in e.g. corre_evaluation_2023
    measure_df['simpson_index'] = 0
    for pathway in pathways:
        measure_df['simpson_index'] = measure_df['simpson_index'] + (
                measure_df[f'mean_identified_as_{pathway}'] * measure_df[f'mean_identified_as_{pathway}'])
    measure_df['simpson_index'] = 1 - measure_df['simpson_index']

    measure_df['number_of_apparent_categories'] = 0
    for pathway in pathways:
        measure_df[f'binary_identified_as_{pathway}'] = 0
        measure_df.loc[measure_df[f'identified_{pathway}_count'] > 0, f'binary_identified_as_{pathway}'] = 1
        measure_df['number_of_apparent_categories'] = measure_df['number_of_apparent_categories'] + measure_df[
            f'binary_identified_as_{pathway}']

    ## Pielou index normalises by the 'richness' for the given genus, i.e. the number of different pathways present
    ## as discussed in corre_evaluation_2023.
    measure_df['pielou_index'] = measure_df['shannon_index'] / (np.log(measure_df['number_of_apparent_categories']))
    measure_df['pielou_index'] = measure_df['pielou_index'].fillna(measure_df['shannon_index'])

    ### More comprehensive normalisation of 'sampling effort' normalises based on the number of tested compounds
    ## This penalises large N.
    measure_df['normalised_simpson'] = measure_df['simpson_index'] * (measure_df['identified_compounds_count']) / (
            measure_df['identified_compounds_count'] - 1)

    measure_df['normalised_shannon'] = measure_df['shannon_index'] / (np.log(measure_df['identified_compounds_count']))

    measure_df.to_csv(os.path.join(processed_compounds_output_path, f'genus_level_pathway_diversity_information.csv'))


def main():
    measure_df = pd.read_csv(os.path.join('outputs','pathways', 'genus_level_pathway_data.csv'), index_col=0)
    get_pathway_based_diversity_measures_for_genera(measure_df, get_pathway_categories())


if __name__ == '__main__':
    main()
