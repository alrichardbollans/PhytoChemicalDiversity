import os
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import spearmanr

from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, RARE_METRICS


def get_sp_working_data():
    trait_data = pd.read_csv(
        os.path.join('..', 'collect_and_compile_data', 'get_diversity_metrics', 'outputs', 'group_data', 'accepted_species_transformed.csv'),
        index_col=0)

    trait_data = trait_data.rename(columns={
        'Assigned_group': 'accepted_species'
    })

    return trait_data


def get_region_counts():
    working_data = get_sp_working_data()
    region_data = pd.read_csv(os.path.join('..', 'collect_and_compile_data/get_diversity_metrics/outputs/group_info/native_regions.csv'), index_col=0)
    species = region_data['accepted_species'].unique()
    out_data = []
    for sp in species:
        species_data = region_data[region_data['accepted_species'] == sp]
        out_data.append([sp, len(species_data['Assigned_group'].unique())])
    out_df = pd.DataFrame(out_data, columns=['accepted_species', 'No. of Regions'])
    out_df = pd.merge(working_data, out_df, how='left', on='accepted_species')
    return out_df


def sp_correlation_calculations(data, metric: str, tag: str, method='spearman'):
    '''
    :param data: scaled data
    :param metric:
    :return:
    '''
    data = data.dropna(subset=['No. of Regions', metric], how='any')
    # Correlation analysis
    # Univariate analyses showing correlations exist
    correlation_matrix = data[['No. of Regions', metric]].corr(method=method)[metric]
    correlation_matrix = correlation_matrix.loc[['No. of Regions']]
    correlation_matrix.columns = [f'{tag}_{metric}']
    print("\nCorrelation Matrix:")
    print(correlation_matrix)

    region_count = spearmanr(data['No. of Regions'], data[metric])

    # just check calculations are same as for plots
    last_computed_value = correlation_matrix.loc['No. of Regions']
    assert round(region_count.correlation, 10) == round(last_computed_value, 10)

    spearmanr_df = pd.DataFrame([region_count.pvalue], index=['No. of Regions'],
                                columns=[f'{tag}_{metric}'])
    return correlation_matrix, spearmanr_df


def plot_distributions(working_data):
    Path(os.path.join('outputs', 'metric_correlations')).mkdir(parents=True, exist_ok=True)
    with sns.plotting_context("notebook", font_scale=2.5):
        ax = sns.pairplot(working_data[METRICS + ['No. of Regions']])
        ax.set(xticklabels=[], yticklabels=[])  # remove the tick labels
        ax.tick_params(bottom=False, left=False)  # remove the ticks
        plt.tight_layout()
        plt.savefig(os.path.join('outputs', 'metric_correlations', 'metric_distributions.png'), dpi=300)
        plt.close()

    # sns.pairplot(working_data[RARE_METRICS + ['No. of Regions']])
    # plt.savefig(os.path.join('outputs', 'metric_correlations', 'rare', 'rare_metric_distributions.png'), dpi=300)
    # plt.close()


def main(metrics, outpath):
    Path(outpath).mkdir(parents=True, exist_ok=True)
    correlations = pd.DataFrame()
    correlations_p = pd.DataFrame()

    working_data = get_region_counts()
    plot_distributions(working_data)
    for metric in metrics:
        tag = 'accepted_species'
        correlation_matrix, spearmanr_df = sp_correlation_calculations(working_data, metric=metric, tag=tag)
        correlations = pd.concat([correlations, correlation_matrix], axis=1)
        correlations_p = pd.concat([correlations_p, spearmanr_df], axis=1)
    correlations.to_csv(os.path.join(outpath, 'correlations.csv'))

    fig, ax = plt.subplots(figsize=(5, 1.5))
    sns.heatmap(correlations.loc[['No. of Regions']], cmap='inferno', annot=True, cbar=True, vmin=0, vmax=1)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, 'correlation_heatmap.jpg'), dpi=300)
    plt.close()

    correlations_p.to_csv(os.path.join(outpath, 'correlations_p.csv'))


if __name__ == '__main__':
    main(METRICS, os.path.join('outputs', 'correlations'))
    # main(RARE_METRICS, os.path.join('outputs', 'correlations', 'rare'))
