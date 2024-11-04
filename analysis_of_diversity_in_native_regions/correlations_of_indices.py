import os
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import spearmanr

from analysis_of_diversity_in_native_regions.analyse_relation_to_pd import get_working_data
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, PATHWAY_INDICES, FAD_INDICES
import seaborn as sns


def plot_distributions():

    sns.pairplot(working_data[METRICS])
    Path(os.path.join('outputs', 'distributions')).mkdir(parents=True, exist_ok=True)
    plt.savefig(os.path.join('outputs', 'metric_correlations', 'metric_distributions.png'), dpi=300)
    plt.close()

    sns.pairplot(working_data[PATHWAY_INDICES + ['GroupSize_Pathways']])
    Path(os.path.join('outputs', 'distributions')).mkdir(parents=True, exist_ok=True)
    plt.savefig(os.path.join('outputs', 'metric_correlations', 'pathway_N_distributions.png'), dpi=300)
    plt.close()


    sns.pairplot(working_data[FAD_INDICES + ['GroupSize_FAD']])
    Path(os.path.join('outputs', 'distributions')).mkdir(parents=True, exist_ok=True)
    plt.savefig(os.path.join('outputs', 'metric_correlations', 'FAD_N_distributions.png'), dpi=300)
    plt.close()

    for h in METRICS:
        sns.pairplot(working_data[[h, 'number_of_species_in_group', 'phylogenetic_diversity']])
        plt.savefig(os.path.join('outputs', 'distributions', h + '_distributions.png'), dpi=300)
        plt.close()


def output_correlations():
    out_dir = os.path.join('outputs', 'metric_correlations')
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    corr_df = working_data[METRICS]
    corr_df = corr_df.corr(method='spearman')
    corr_df.to_csv(os.path.join(out_dir, 'correlations.csv'))
    # plot the heatmap

    plot_corr_df = corr_df.drop(columns=['APWD'])
    plot_corr_df = plot_corr_df.drop('H')
    # Getting the Upper Triangle of the co-relation matrix
    mask = np.zeros_like(plot_corr_df)
    mask[np.triu_indices_from(mask)] = True

    # Want diagonal elements as well
    mask[np.diag_indices_from(mask)] = False
    sns.heatmap(plot_corr_df, cmap='coolwarm', annot=True, mask=mask, cbar=False)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'diversity_heatmap.jpg'), dpi=300)
    plt.close()

    dfcols = working_data[METRICS]
    pvalues = pd.DataFrame(index=METRICS, columns=METRICS)
    for r in METRICS:
        for c in METRICS:
            w = dfcols.dropna(subset=[r, c],how='any')
            sp = spearmanr(w[r], w[c])
            # just check calculations are same as for plots
            last_computed_value = corr_df.loc[r][c]
            assert round(sp.correlation,10) == round(last_computed_value,10)
            pvalues[r][c] = sp.pvalue

    pvalues.to_csv(os.path.join(out_dir, 'correlations_pvalues.csv'))



def main():
    # plot_distributions()
    output_correlations()

if __name__ == '__main__':
    working_data = get_working_data()

    main()
