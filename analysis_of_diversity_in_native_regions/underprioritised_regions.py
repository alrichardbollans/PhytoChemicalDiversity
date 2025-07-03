import os.path
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from loess.loess_2d import loess_2d
from matplotlib import pyplot as plt
import statsmodels.api as sm
from sklearn.metrics import r2_score
from sklearn.preprocessing import PolynomialFeatures
from analysis_of_diversity_in_native_regions.helper_functions import get_working_data
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, RARE_METRICS

outdir = os.path.join('outputs', 'underprioritised_regions')
os.makedirs(outdir, exist_ok=True)


def get_reg_data(metric):
    reg_data = working_data[['Group', 'PD', 'SR', metric]].dropna(
        subset=['PD', 'SR', metric], how='any')

    return reg_data


def multiple_loess(metric):
    # From Cappellari et al., ‘The ATLAS3D Project – XX. Mass–Size and Mass–σ Distributions of Early-Type Galaxies’.
    # 2D version is based on Cleveland and Devlin, ‘Locally Weighted Regression: An Approach to Regression Analysis by Local Fitting’.
    reg_data = get_reg_data(metric)
    x, y, z = reg_data['PD'].values, reg_data['SR'].values, reg_data[metric].values
    zout, wout = loess_2d(x, y, z, xnew=None, ynew=None, degree=1, frac=0.5,
                          npoints=None, rescale=False, sigz=None)
    r2 = r2_score(z, zout)

    analysis_df = reg_data
    analysis_df['expected_diversity'] = zout
    analysis_df['wout'] = wout
    analysis_df['r2'] = r2
    analysis_df[f'residuals'] = analysis_df[metric] - analysis_df['expected_diversity']

    std_residual = analysis_df[f'residuals'].std()
    mean_residual = analysis_df[f'residuals'].mean()
    analysis_df['highlight_high'] = analysis_df[f'residuals'] > ((2 * std_residual) + mean_residual)
    analysis_df['highlight_low'] = analysis_df[f'residuals'] < (mean_residual - (2 * std_residual))

    return zout, wout, r2, analysis_df


def plot_2d_annotated_regression_data(data, metric_outpath, x_var, y_var):
    # Set up the plot
    import seaborn as sns
    from pypalettes import load_cmap
    # plt.figure(figsize=(10, 6))
    # sns.set_style("whitegrid")
    # Scatter plot for APWD vs phylogenetic_diversity
    # colors = load_cmap('inferno').hex
    # print(colors)
    data['color'] = np.where((data['highlight_high'] == True), '#d12020', 'grey')
    data['color'] = np.where((data['highlight_low'] == True), '#5920ff', data['color'])

    sns.scatterplot(x=x_var, y=y_var, data=data, color=data['color'], edgecolor="black", alpha=0.8)

    # Highlight points where 'highlight' is True

    highlighted_data = data[(data['highlight_high'] == True) | (data['highlight_low'] == True)]

    for _, row in highlighted_data.iterrows():
        upshift = 0
        left_shift = 0
        # if row['highlight_high']:  # in ['COR', 'KZN', 'ALD', 'AZO', 'ROD', 'SEY', 'LDV']:
        #     upshift = 0
        #     left_shift = 0.05
        #     if row['Group'] in ['ALD', 'SEY']:
        #         upshift = 0.09
        #     if row['Group'] in ['AZO']:
        #         upshift = -0.2
        #     if row['Group'] in ['SEY']:
        #         left_shift = 0
        #
        #     if row['Group'] in ['COR']:
        #         left_shift = -0.05
        #         upshift = 0.11
        # if row['highlight_low']:
        #     upshift = 0
        #     left_shift = -0.45
        #     if row['Group'] in ['CLS']:
        #         left_shift = -0.38
        #     if row['Group'] in ['CLS', 'MSO']:
        #         upshift = -0.2
        #     if row['Group'] in ['AGS']:
        #         upshift = -0.13
        #     if row['Group'] in ['NZN']:
        #         upshift = 0.04
        #         left_shift = -0.50

        plt.annotate(row['Group'], (row[x_var] + left_shift, row[y_var] + upshift), ha='right', color='black')

    # Line plot for expected_diversity vs xvar

    ## Add estimator to smooth cases where multiple observations of the y variable at the same x
    sns.lineplot(x=x_var, y='expected_diversity', estimator='mean', color='black', linestyle='--', data=data)

    # Labels and legend

    plt.xlabel(x_var)

    plt.ylabel(y_var)

    plt.title('')

    # plt.legend(loc='upper right')

    plt.savefig(os.path.join(metric_outpath, f'{x_var}_plot_with_outliers.jpg'), dpi=300)
    plt.close()
    plt.clf()
    sns.reset_orig()


def find_regression_model(metric):
    regression_out_dir = os.path.join(outdir, 'regressions', metric)
    os.makedirs(regression_out_dir, exist_ok=True)
    reg_data = get_reg_data(metric)
    y_var = metric
    x_var = 'PD'
    # Step 1: Fit a linear regression model
    X = reg_data[[x_var]].values  # Independent variable
    X_to_plot = reg_data[x_var].values  # Independent variable (SR)
    # scaled_data[metric] = np.log(scaled_data[metric])
    y = reg_data[y_var].values  # Dependent variable (diversity)

    data = []

    # Fit LOESS with outlier robustness (iterations downweight outliers)
    loess_prediction = sm.nonparametric.lowess(exog=X_to_plot, endog=y, return_sorted=False)
    expected_diversity = loess_prediction
    residuals = y - expected_diversity
    # Calculate R² (coefficient of determination)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)
    print(f"R² (Coefficient of Determination): {r_squared:.3f}")

    sns.scatterplot(x=X_to_plot, y=y, edgecolor="black", alpha=0.8)
    sns.lineplot(x=X_to_plot, y=expected_diversity, color='black', linestyle='--')
    plt.xlabel('Phylogenetic Diversity')
    plt.ylabel(metric)
    plt.savefig(os.path.join(regression_out_dir, 'LOESS.jpg'), dpi=300)
    plt.close()
    data.append([f'LOESS', r_squared])

    df = pd.DataFrame(data, columns=['Model', 'R-squared'])
    df.to_csv(os.path.join(regression_out_dir, 'model_comparison.csv'))
    return loess_prediction


def main():
    collected_highlights = {}
    for m in METRICS + RARE_METRICS:
        metric_outpath = os.path.join(outdir, 'regressions', m)
        Path(metric_outpath).mkdir(parents=True, exist_ok=True)

        zout, wout, r2, analysis_df = multiple_loess(m)
        analysis_df.to_csv(os.path.join(metric_outpath, f'loess_outputs.csv'))
        plot_2d_annotated_regression_data(analysis_df, metric_outpath, 'PD', m)
        plot_2d_annotated_regression_data(analysis_df, metric_outpath, 'SR', m)

        highlights = analysis_df[(analysis_df['highlight_high'] == True)]['Group'].tolist()
        collected_highlights[m] = highlights
    consistent_compound_highlights = set(collected_highlights['APWD'])
    for m in ['APWD']:
        consistent_compound_highlights = consistent_compound_highlights.intersection(collected_highlights[m])

    print(consistent_compound_highlights)
    consistent_pathway_highlights = set(collected_highlights['H'])
    for m in ['H', 'Hbc', 'G']:
        consistent_pathway_highlights = consistent_pathway_highlights.intersection(collected_highlights[m])
    print(consistent_pathway_highlights)


if __name__ == '__main__':
    working_data = get_working_data()

    main()
