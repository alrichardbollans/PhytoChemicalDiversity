import os

import pandas as pd
from scipy.stats import spearmanr

from helper_functions import get_working_data_for_medicinal_species_analysis
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, RARE_METRICS


def duke():
    PD_indep_data = pd.read_csv(os.path.join('outputs', 'PD_SR_regression', 'SR-Independent PD.csv'))[['Group', "PD'"]]

    working_data = get_working_data_for_medicinal_species_analysis()
    test_len = len(working_data)
    working_data = pd.merge(working_data, PD_indep_data, how='left', on='Group')

    medicinal_data = pd.read_csv(
        os.path.join('..','..', 'collect_and_compile_data', 'compile_trait_data', 'outputs', 'medicinal_activities_per_region.csv'))[
        ['Assigned_group', 'number_medicinal_activities_in_region']].rename(columns={'Assigned_group': 'Group'})
    data = pd.merge(working_data, medicinal_data, on='Group')
    data = data.dropna(subset=['PD', 'number_medicinal_activities_in_region'], how='any')

    for m in ['PD', "PD'", 'SR']+METRICS + RARE_METRICS:
        apwd = spearmanr(data[m], data['number_medicinal_activities_in_region'])
        print(f'{m}: {apwd.correlation}')
    pass
if __name__ == '__main__':
    duke()