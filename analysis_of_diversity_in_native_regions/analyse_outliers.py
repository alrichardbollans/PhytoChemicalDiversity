import os.path
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import stats
from sklearn.linear_model import LinearRegression
from analysis_of_diversity_in_native_regions.analyse_relation_to_pd import get_working_data
from analysis_of_diversity_in_native_regions.spatial_plots import plot_dist_of_metric
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, RARE_METRICS


# def outliers():
#     working_data = get_working_data()[['Group'] + METRICS]
#     Path(os.path.join('outputs', 'outliers')).mkdir(parents=True, exist_ok=True)
#
#     for metric in METRICS:
#         metric_data = working_data[['Group', metric]]
#         # IQR Method (Interquartile Range)
#         # This method doesn't assume normality at all.
#         # It looks at the middle 50% of the data (between the 1st and 3rd quartile) and identifies outliers based on how far a value falls outside this range.
#         # # Calculate Q1 (25th percentile) and Q3 (75th percentile)
#         # Q1 = metric_data[metric].quantile(0.25)
#         # Q3 = metric_data[metric].quantile(0.75)
#         # IQR = Q3 - Q1
#         #
#         # # Define outliers using the IQR method
#         # lower_bound = Q1 - 1.5 * IQR
#         # upper_bound = Q3 + 1.5 * IQR
#
#         lower_bound = metric_data[metric].quantile(0.05)
#         upper_bound = metric_data[metric].quantile(0.95)
#
#         positive_outliers = metric_data[(metric_data[metric] > upper_bound)]
#         neg_outliers = metric_data[(metric_data[metric] < lower_bound)]
#         positive_outliers.to_csv(os.path.join('outputs', 'outliers', f'{metric}_high.csv'))
#         neg_outliers.to_csv(os.path.join('outputs', 'outliers', f'{metric}_low.csv'))


def using_models():
    for metric in METRICS + RARE_METRICS:
        if metric in METRICS:
            outpath = os.path.join('outputs', 'regressions')
        else:
            outpath = os.path.join('outputs', 'regressions', 'rare')
        Path(outpath).mkdir(parents=True, exist_ok=True)
        working_data = get_working_data()[['Group', 'Phylogenetic Diversity', 'number_of_species_in_group', metric]]
        working_data = working_data.dropna(subset=['Phylogenetic Diversity', metric])
        # Step 1: Fit a linear regression model
        X = working_data[['Phylogenetic Diversity']].values  # Independent variable (species richness)
        # scaled_data[metric] = np.log(scaled_data[metric])
        y = working_data[metric].values  # Dependent variable (diversity)

        model = LinearRegression()
        model.fit(X, y)

        r_squared = model.score(X, y)
        print(f'R-squared value: {r_squared:.4f}')
        working_data['r_squared'] = r_squared
        # Step 2: Predict the expected diversity based on species richness
        working_data['expected_diversity'] = model.predict(X)

        # Step 3: Calculate the residuals (observed - expected)
        working_data[f'{metric}_residuals'] = working_data[metric] - working_data['expected_diversity']

        # Step 4: Highlight cases with large residuals
        # Let's consider residuals greater than 2 standard deviations as "large differences"

        std_residual = working_data[f'{metric}_residuals'].std()
        mean_residual = working_data[f'{metric}_residuals'].mean()
        if mean_residual > 0.0000001 or mean_residual < -0.0000001:
            raise ValueError
        working_data['highlight_high'] = working_data[f'{metric}_residuals'] > (2 * std_residual)
        working_data['highlight_low'] = working_data[f'{metric}_residuals'] < -(2 * std_residual)

        # Print the cases with large differences
        # print("Cases with large differences (residuals):")
        # print(working_data[working_data['highlight']])

        working_data.to_csv(os.path.join(outpath, f'{metric}.csv'))

        def estimator(x):
            return model.predict(np.array(x).reshape(-1, 1))

        plot_annotated_regression_data(working_data, outpath, metric, 'Phylogenetic Diversity')

        plot_dist_of_metric(working_data, f'{metric}_residuals', out_path=os.path.join(outpath, 'dists', f'{metric}.jpg'))


def plot_annotated_regression_data(data, outpath, metric, x_var):
    """

    Plots a scatter plot of 'phylogenetic_diversity' vs 'APWD' and a line plot of

    'phylogenetic_diversity' vs 'expected_diversity'. Highlights data points with labels

    where 'highlight' is True.

    """

    # Set up the plot
    import seaborn as sns
    from pypalettes import load_cmap
    # plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")
    # Scatter plot for APWD vs phylogenetic_diversity
    colors = load_cmap('Acadia').hex
    data['color'] = np.where((data['highlight_high'] == True) | (data['highlight_low'] == True), colors[0], colors[2])

    sns.scatterplot(x=x_var, y=metric, data=data, color=data['color'], edgecolor="black", alpha=0.8)

    # Highlight points where 'highlight' is True

    highlighted_data = data[(data['highlight_high'] == True) | (data['highlight_low'] == True)]

    for _, row in highlighted_data.iterrows():
        if row['Group'] in ['COR', 'KZN', 'ALD', 'AZO', 'ROD', 'SEY', 'LDV']:
            upshift =0
            left_shift = 0.05
            if row['Group'] in ['ALD', 'SEY']:
                upshift = 0.09
            if row['Group'] in ['AZO']:
                upshift = -0.2
            if row['Group'] in [ 'SEY']:
                left_shift = 0

            if row['Group'] in [ 'COR']:
                left_shift = -0.05
                upshift = 0.11

            plt.annotate(row['Group'], (row[x_var]-left_shift, row[metric]+upshift), ha='right', color='black')

    # Line plot for expected_diversity vs phylogenetic_diversity

    sns.lineplot(x=x_var, y='expected_diversity', color='black', linestyle='--', data=data)

    # Labels and legend

    plt.xlabel(x_var)

    plt.ylabel(metric)

    plt.title('')

    # plt.legend(loc='upper right')

    plt.savefig(os.path.join(outpath, f'{metric}.jpg'), dpi=300)
    plt.close()
    sns.reset_orig()


def collect_highlights():
    out_dict = {}
    fileNames = os.listdir(os.path.join('outputs', 'regressions'))
    rare_fileNames = os.listdir(os.path.join('outputs', 'regressions', 'rare'))
    for f in fileNames + rare_fileNames:
        if f.endswith(".csv"):
            if f in fileNames:
                df = pd.read_csv(os.path.join('outputs', 'regressions', f))
            else:
                df = pd.read_csv(os.path.join('outputs', 'regressions', 'rare', f))

            highlights = df[df['highlight_high'] == True]['Group'].tolist()
            lowlights = df[df['highlight_low'] == True]['Group'].tolist()
            out_dict[f[:-4]] = {'high': highlights, 'low': lowlights}

    interesting_highlights = set(out_dict['APWD']['high'])
    for m in ['APWD']:
        interesting_highlights = interesting_highlights.intersection(set(out_dict[m]['high']))
    print(f'Interesting highlights: {interesting_highlights}')
    universal_highlights = set(out_dict[METRICS[0]]['high'])
    for m in METRICS:
        universal_highlights = universal_highlights.intersection(set(out_dict[m]['high']))

    universal_lowlights = set(out_dict[METRICS[0]]['low'])
    for m in METRICS:
        universal_lowlights = universal_lowlights.intersection(set(out_dict[m]['low']))

    print(universal_highlights)
    print(universal_lowlights)

    universal_highlights_without_J = set(out_dict['APWD']['high'])
    for m in METRICS:
        if m != 'J':
            universal_highlights_without_J = universal_highlights_without_J.intersection(set(out_dict[m]['high']))

    universal_lowlights_without_J = set(out_dict['APWD']['low'])
    for m in METRICS:
        if m != 'J':
            universal_lowlights_without_J = universal_lowlights_without_J.intersection(set(out_dict[m]['low']))

    print(universal_highlights)
    print(universal_lowlights)
    return out_dict


def get_info_about_a_region(region_code: str):
    # Number of species
    # sps names

    sp_data = pd.read_csv(os.path.join('..', 'collect_and_compile_data/get_diversity_metrics/outputs/group_info/native_regions.csv'), index_col=0)
    region_data = sp_data[sp_data['Assigned_group'] == region_code]
    species = region_data['accepted_species'].unique()
    num_species = len(species)

    endemic_species = []
    for sp in species:
        if len(sp_data[sp_data['accepted_species'] == sp]['Assigned_group'].unique().tolist()) == 1:
            endemic_species.append(sp)

    # Number of compounds
    group_data = pd.read_csv(os.path.join('..', 'collect_and_compile_data/get_diversity_metrics/outputs/group_data/native_regions.csv'), index_col=0)
    group_data = group_data[group_data['Assigned_group'] == region_code]
    assert len(group_data) == 1
    num_species1 = group_data['number_of_species_in_group'].iloc[0]
    assert num_species == num_species1
    num_compounds = group_data['GroupSize_FAD'].iloc[0]

    # highlights for which measures

    highlight_data = collect_highlights()
    highlight_metrics = []
    for m in highlight_data:
        if region_code in highlight_data[m]['high']:
            highlight_metrics.append(m)

    # Maybe a tree plot would be nice?
    # Phydiv quartile.
    all_phydiv_df = pd.read_csv(os.path.join('..', 'collect_and_compile_data/get_phylogenetic_diversities/outputs/group_data/native_regions.csv'))
    describe = all_phydiv_df.describe()
    phydiv_df = all_phydiv_df[all_phydiv_df['Group'] == region_code]
    assert len(phydiv_df) == 1
    phydiv = phydiv_df['phylogenetic_diversity'].iloc[0]
    percentile = stats.percentileofscore(all_phydiv_df['phylogenetic_diversity'].dropna().values, phydiv)

    formatted_species = [x.replace(' ', '_') for x in species]

    out_df = pd.DataFrame([num_species, str(species),str(endemic_species), formatted_species,num_compounds, str(highlight_metrics), phydiv, percentile],
                          index=['num_species', 'species','endemic_species','formatted_species', 'num_compounds', 'highlight_metrics', 'phydiv', 'percentile'], columns=[region_code])
    out_df.to_csv(os.path.join('outputs', 'outlier_info', f'{region_code}.csv'))


if __name__ == '__main__':
    using_models()
    collect_highlights()
    get_info_about_a_region('ALD') # Aldabra
    get_info_about_a_region('ROD') # Rodrigues
    get_info_about_a_region('COR') # Corsica
    get_info_about_a_region('LDV') # Laccadive Is.
    get_info_about_a_region('KZN') # Kazan-retto
    get_info_about_a_region('AZO') # Azores
    get_info_about_a_region('SEY') # Seychelles
