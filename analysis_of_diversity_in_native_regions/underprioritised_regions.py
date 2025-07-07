import os.path
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from matplotlib import pyplot as plt

from analysis_of_diversity_in_native_regions.helper_functions import get_working_data
from collect_and_compile_data.collect_compound_data import COMPOUND_ID_COL
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, RARE_METRICS

outdir = os.path.join('outputs', 'underprioritised_regions')
os.makedirs(outdir, exist_ok=True)


def get_reg_data(metric):
    reg_data = working_data[['Group', 'PD', 'SR', metric]].dropna(
        subset=['PD', 'SR', metric], how='any')

    return reg_data


def plot_2d_annotated_regression_data(data, metric_outpath, x_var, y_var, extras_to_annotate:list=None):
    # Set up the plot
    import seaborn as sns
    # from pypalettes import load_cmap
    # plt.figure(figsize=(10, 6))
    # sns.set_style("whitegrid")
    # Scatter plot for APWD vs phylogenetic_diversity
    # colors = load_cmap('inferno').hex
    # print(colors)
    data['color'] = np.where((data[f'{y_var}_highlight_high'] == True), '#d12020', 'grey')
    data['color'] = np.where((data[f'{y_var}_highlight_low'] == True), '#5920ff', data['color'])

    sns.scatterplot(x=x_var, y=y_var, data=data, color=data['color'], edgecolor="black", alpha=0.8)

    # Highlight points where 'highlight' is True

    highlighted_data = data[(data[f'{y_var}_highlight_high'] == True)|(data[f'{y_var}_highlight_low'] == True)]
    to_annotate = highlighted_data[f'Group'].unique().tolist()

    if extras_to_annotate is not None:
        to_annotate += extras_to_annotate
    for _, row in data.iterrows():
        if row['Group'] in to_annotate:
            upshift = 0
            left_shift = -0.05
            if row['Group'] in ['FIJ']:
                left_shift = 0.3
            if row['Group'] in ['VAN', 'NGR']:
                left_shift = 0.45
            if row['Group'] in ['TON', 'TCI', 'KAN']:
                upshift = -0.18
            if row['Group'] == 'SAM':
                upshift = 0.05
            plt.annotate(row['Group'], (row[x_var] + left_shift, row[y_var] + upshift), ha='right', color='black')

    # Line plot for expected_diversity vs xvar

    ## Add estimator to smooth cases where multiple observations of the y variable at the same x
    sns.lineplot(x=x_var, y=f'{y_var}_loess_prediction', estimator='mean', color='black', linestyle='--', data=data)

    # Labels and legend

    plt.xlabel('Phylogenetic Diversity')

    plt.ylabel(y_var, rotation=0)

    plt.title('')

    # plt.legend(loc='upper right')

    plt.savefig(os.path.join(metric_outpath, f'{x_var}_plot_with_outliers.jpg'), dpi=300)
    plt.close()
    plt.clf()
    sns.reset_orig()


def get_loess_outputs(metric):
    regression_out_dir = os.path.join(outdir, 'regressions', metric)
    os.makedirs(regression_out_dir, exist_ok=True)
    reg_data = get_reg_data(metric)
    y_var = metric
    x_var = 'PD'
    # X = reg_data[[x_var]].values  # Independent variable
    X_to_plot = reg_data[x_var].values  # Independent variable
    # scaled_data[metric] = np.log(scaled_data[metric])
    y = reg_data[y_var].values  # Dependent variable

    # Fit LOESS with outlier robustness (iterations downweight outliers)
    loess_prediction = sm.nonparametric.lowess(exog=X_to_plot, endog=y, return_sorted=False)
    expected_diversity = loess_prediction
    residuals = y - expected_diversity
    # Calculate R² (coefficient of determination)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f"R² (Coefficient of Determination): {r_squared:.3f}")

    # sns.scatterplot(x=X_to_plot, y=y, edgecolor="black", alpha=0.8)
    # sns.lineplot(x=X_to_plot, y=expected_diversity, color='black', linestyle='--')
    # plt.xlabel('Phylogenetic Diversity')
    # plt.ylabel(metric)
    # plt.savefig(os.path.join(regression_out_dir, 'LOESS.jpg'), dpi=300)
    # plt.close()

    out_df = reg_data
    out_df[f'{metric}_loess_prediction'] = loess_prediction
    out_df[f'{metric}_residuals'] = out_df[metric] - out_df[f'{metric}_loess_prediction']

    std_residual = out_df[f'{metric}_residuals'].std()
    mean_residual = out_df[f'{metric}_residuals'].mean()
    out_df[f'{metric}_highlight_high'] = out_df[f'{metric}_residuals'] > ((2 * std_residual) + mean_residual)
    out_df[f'{metric}_highlight_low'] = out_df[f'{metric}_residuals'] < (mean_residual - (2 * std_residual))
    out_df['R2'] = r_squared
    out_df.to_csv(os.path.join(regression_out_dir, 'loess_outputs.csv'))
    extras_to_annotate = None
    if metric == 'J_Rare':

        extras_to_annotate = ['AGW', 'FIJ', 'NUE', 'SAM', 'TAS', 'TON', 'VAN']
    # elif metric == 'J':
    #     extras_to_annotate = ['NCS', 'TZK', 'UZB']
    plot_2d_annotated_regression_data(out_df, regression_out_dir, 'PD', metric, extras_to_annotate=extras_to_annotate)

    return out_df


def get_info_about_a_region(region, metric):
    out_data = []
    all_region_data = pd.read_csv(
        os.path.join('..', 'collect_and_compile_data', 'get_diversity_metrics', 'outputs', 'group_info',
                     'native_regions.csv'))
    region_data = all_region_data[all_region_data['Assigned_group'] == region]

    species_in_region = region_data['accepted_species'].unique().tolist()
    all_compound_data = pd.read_csv(
        os.path.join('..', 'collect_and_compile_data', 'collect_compound_data', 'outputs',
                     'all_species_compound_data.csv'))[['accepted_species', COMPOUND_ID_COL]]
    all_compound_data = all_compound_data.drop_duplicates(subset=['accepted_species', COMPOUND_ID_COL])

    region_compound_data = all_compound_data[all_compound_data['accepted_species'].isin(species_in_region)]
    number_compounds_from_region = len(region_compound_data[COMPOUND_ID_COL].unique().tolist())

    out_data.append([number_compounds_from_region])
    out_data.append([str(species_in_region)])
    out_data.append([len(species_in_region)])
    extra_cols = []
    for sp in species_in_region:
        species_compounds = region_compound_data[region_compound_data['accepted_species'] == sp]
        number_compounds_from_sp = len(species_compounds[COMPOUND_ID_COL].unique().tolist())
        out_data.append([number_compounds_from_sp])
        extra_cols.append(f'number_compounds_from_{sp}')
        species_regions = all_region_data[all_region_data['accepted_species'] == sp]
        number_regions_for_species = len(species_regions['Assigned_group'].unique().tolist())
        out_data.append([number_regions_for_species])
        extra_cols.append(f'number_regions_for_{sp}')

    out_df = pd.DataFrame(out_data, index=['number_compounds_from_region', 'species_in_region', 'number species_in_region'] + extra_cols)
    outpath = os.path.join(outdir, 'region_info', metric)
    Path(outpath).mkdir(parents=True, exist_ok=True)
    out_df.to_csv(os.path.join(outpath, f'{region}_info.csv'))


def main():
    collected_highlights = {}
    for m in METRICS + RARE_METRICS:
        metric_outpath = os.path.join(outdir, 'regressions', m)
        Path(metric_outpath).mkdir(parents=True, exist_ok=True)

        metric_analysis_df = get_loess_outputs(m)

        highlights = metric_analysis_df[(metric_analysis_df[f'{m}_highlight_high'] == True)]['Group'].tolist()
        collected_highlights[m] = highlights
        for r in highlights + metric_analysis_df[(metric_analysis_df[f'{m}_highlight_low'] == True)]['Group'].tolist():
            get_info_about_a_region(r, m)
    print(collected_highlights)
    consistent_compound_highlights = set(collected_highlights['APWD'])
    for m in ['FAD', 'MFAD', 'APWD']:
        consistent_compound_highlights = consistent_compound_highlights.intersection(collected_highlights[m])
    print(consistent_compound_highlights)

    consistent_rare_compound_highlights = set(collected_highlights['APWD_Rare'])
    for m in ['FAD_Rare', 'MFAD_Rare', 'APWD_Rare']:
        consistent_rare_compound_highlights = consistent_rare_compound_highlights.intersection(collected_highlights[m])
    print(consistent_rare_compound_highlights)

    consistent_pathway_highlights = set(collected_highlights['H'])
    for m in ['H', 'Hbc', 'G']:
        consistent_pathway_highlights = consistent_pathway_highlights.intersection(collected_highlights[m])
    print(consistent_pathway_highlights)

    consistent_rare_pathway_highlights = set(collected_highlights['G_Rare'])
    for m in ['H_Rare', 'Hbc_Rare', 'G_Rare']:
        consistent_rare_pathway_highlights = consistent_rare_pathway_highlights.intersection(collected_highlights[m])
    print(consistent_rare_pathway_highlights)

    out_df = pd.DataFrame([str(consistent_compound_highlights), str(consistent_rare_compound_highlights),
                           str(consistent_pathway_highlights), str(consistent_rare_pathway_highlights), str(collected_highlights['J']), str(collected_highlights['J_Rare'])],
                          index=['compounds', 'rare_compounds', 'pathways', 'rare_pathways', 'J', 'rare_J'])
    out_df.to_csv(os.path.join(outdir, 'highlights.csv'))


if __name__ == '__main__':
    working_data = get_working_data()

    main()
