import os
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import spearmanr

from analysis_of_diversity_in_native_regions.analyse_relation_to_pd import get_working_data
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, PATHWAY_INDICES, FAD_INDICES, RARE_METRICS
import seaborn as sns


def plot_distributions():
    Path(os.path.join('outputs', 'metric_correlations')).mkdir(parents=True, exist_ok=True)
    sns.pairplot(working_data[METRICS + ['N']])
    plt.savefig(os.path.join('outputs', 'metric_correlations', 'metric_distributions.png'), dpi=300)
    plt.close()

    sns.pairplot(working_data[RARE_METRICS + ['N', 'GroupSize_Pathways']])
    plt.savefig(os.path.join('outputs', 'metric_correlations', 'rare', 'rare_metric_distributions.png'), dpi=300)
    plt.close()

    sns.pairplot(working_data[PATHWAY_INDICES + ['GroupSize_Pathways']])
    plt.savefig(os.path.join('outputs', 'metric_correlations', 'pathway_N_distributions.png'), dpi=300)
    plt.close()


    sns.pairplot(working_data[FAD_INDICES + ['N']])
    plt.savefig(os.path.join('outputs', 'metric_correlations', 'FAD_N_distributions.png'), dpi=300)
    plt.close()

    for h in METRICS:
        sns.pairplot(working_data[[h, f'{h}_Rare', 'number_of_species_in_group', 'Phylogenetic Diversity', 'N']])
        plt.savefig(os.path.join('outputs', 'distributions', h + '_distributions.png'), dpi=300)
        plt.close()


def output_correlations(metrics, out_dir):
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    corr_df = working_data[metrics]
    corr_df = corr_df.corr(method='spearman')
    corr_df.to_csv(os.path.join(out_dir, 'correlations.csv'))
    # plot the heatmap

    plot_corr_df = corr_df.drop(columns=[corr_df.columns[-1]])
    plot_corr_df = plot_corr_df.drop(corr_df.columns[0])
    # Getting the Upper Triangle of the co-relation matrix
    mask = np.zeros_like(plot_corr_df)
    mask[np.triu_indices_from(mask)] = True

    # Want diagonal elements as well
    mask[np.diag_indices_from(mask)] = False
    sns.heatmap(plot_corr_df, cmap='viridis', annot=True, mask=mask, cbar=False)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'diversity_heatmap.jpg'), dpi=300)
    plt.close()

    dfcols = working_data[metrics]
    pvalues = pd.DataFrame(index=metrics, columns=metrics)
    for r in metrics:
        for c in metrics:
            w = dfcols.dropna(subset=[r, c],how='any')
            sp = spearmanr(w[r], w[c])
            # just check calculations are same as for plots
            last_computed_value = corr_df.loc[r][c]
            assert round(sp.correlation,10) == round(last_computed_value,10)
            pvalues[r][c] = sp.pvalue

    pvalues.to_csv(os.path.join(out_dir, 'correlations_pvalues.csv'))



def main():
    plot_distributions()
    output_correlations(METRICS, os.path.join('outputs', 'metric_correlations'))
    output_correlations(RARE_METRICS, os.path.join('outputs', 'metric_correlations', 'rare'))

if __name__ == '__main__':
    working_data = get_working_data()
    working_data = working_data.rename(columns={'GroupSize_FAD': 'N'})

    main()
