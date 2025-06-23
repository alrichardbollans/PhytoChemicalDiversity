import os

import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

from analysis_of_diversity_in_native_regions.helper_functions import get_working_data
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, RARE_METRICS


def correlation_calculations(data, metric: str, tag: str, method='spearman'):
    '''
    :param data: scaled data
    :param metric:
    :return:
    '''
    data = data.dropna(subset=['PD','SR', "PD'", metric], how='any')
    # Correlation analysis
    # Univariate analyses showing correlations exist
    correlation_matrix = data[['PD', 'SR',"PD'", metric]].corr(method=method)[metric]
    correlation_matrix = correlation_matrix.loc[['PD', 'SR',"PD'"]]
    correlation_matrix.columns = [f'{tag}_{metric}']
    print("\nCorrelation Matrix:")
    print(correlation_matrix)

    phydiv = spearmanr(data['PD'], data[metric])
    SR = spearmanr(data['SR'], data[metric])
    PD_indep = spearmanr(data["PD'"], data[metric])

    # just check calculations are same as for plots
    last_computed_value = correlation_matrix.loc['PD']
    assert round(phydiv.correlation, 10) == round(last_computed_value, 10)

    spearmanr_df = pd.DataFrame([phydiv.pvalue, SR.pvalue, PD_indep.pvalue], index=['PD', 'SR', "PD'"],
                                columns=[f'{tag}_{metric}'])
    return correlation_matrix, spearmanr_df


def main(metrics, outpath):
    Path(outpath).mkdir(parents=True, exist_ok=True)
    correlations = pd.DataFrame()
    correlations_p = pd.DataFrame()
    PD_indep_data = pd.read_csv(os.path.join('outputs','PD_SR_regression','SR-Independent PD.csv'))[['Group', "PD'"]]
    for metric in metrics:
        tag = 'native_regions'
        working_data = get_working_data()
        test_len = len(working_data)
        working_data = pd.merge(working_data, PD_indep_data, how='left', on='Group')
        assert len(working_data) == test_len
        # scaled_data = get_group_data(metric)
        correlation_matrix, spearmanr_df = correlation_calculations(working_data, metric=metric, tag=tag)
        correlations = pd.concat([correlations, correlation_matrix], axis=1)
        correlations_p = pd.concat([correlations_p, spearmanr_df], axis=1)
    correlations.to_csv(os.path.join(outpath, 'correlations.csv'))

    fig, ax = plt.subplots(figsize=(5, 1.5))
    sns.heatmap(correlations.loc[['PD', 'SR', "PD'"]],cmap='inferno', annot=True, cbar=True, vmin=0, vmax=1)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, 'correlation_heatmap.jpg'), dpi=300)
    plt.close()

    correlations_p.to_csv(os.path.join(outpath, 'correlations_p.csv'))


if __name__ == '__main__':
    main(METRICS, os.path.join('outputs', 'correlations'))
    main(RARE_METRICS, os.path.join('outputs', 'correlations', 'rare'))