import os.path
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression

from analysis_of_all_traits.analyse_groups import get_working_data
from analysis_of_all_traits.spatial_plots import plot_dist_of_metric
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS


def plot_distributions():
    import seaborn as sns
    working_data = get_working_data()
    sns.pairplot(working_data[METRICS])
    Path(os.path.join('outputs', 'distributions')).mkdir(parents=True, exist_ok=True)
    plt.savefig(os.path.join('outputs', 'distributions', 'metric_distributions.png'), dpi=300)
    plt.close()

    for h in METRICS:
        sns.pairplot(working_data[[h, 'number_of_species_in_group', 'phylogenetic_diversity']])
        plt.savefig(os.path.join('outputs', 'distributions', h + '_distributions.png'), dpi=300)
        plt.close()


def outliers():
    working_data = get_working_data()[['Group'] + METRICS]
    Path(os.path.join('outputs', 'outliers')).mkdir(parents=True, exist_ok=True)

    for metric in METRICS:
        metric_data = working_data[['Group', metric]]
        # IQR Method (Interquartile Range)
        # This method doesn't assume normality at all.
        # It looks at the middle 50% of the data (between the 1st and 3rd quartile) and identifies outliers based on how far a value falls outside this range.
        # # Calculate Q1 (25th percentile) and Q3 (75th percentile)
        # Q1 = metric_data[metric].quantile(0.25)
        # Q3 = metric_data[metric].quantile(0.75)
        # IQR = Q3 - Q1
        #
        # # Define outliers using the IQR method
        # lower_bound = Q1 - 1.5 * IQR
        # upper_bound = Q3 + 1.5 * IQR

        lower_bound = metric_data[metric].quantile(0.05)
        upper_bound = metric_data[metric].quantile(0.95)

        positive_outliers = metric_data[(metric_data[metric] > upper_bound)]
        neg_outliers = metric_data[(metric_data[metric] < lower_bound)]
        positive_outliers.to_csv(os.path.join('outputs', 'outliers', f'{metric}_high.csv'))
        neg_outliers.to_csv(os.path.join('outputs', 'outliers', f'{metric}_low.csv'))


def using_models():
    x_vars = ['number_of_species_in_group', 'phylogenetic_diversity']
    for x_var in x_vars:
        outpath = os.path.join('outputs', 'regressions', x_var)
        Path(outpath).mkdir(parents=True, exist_ok=True)

        for metric in METRICS:
            working_data = get_working_data()[['Group', 'phylogenetic_diversity', 'number_of_species_in_group', 'N', metric]]
            working_data = working_data.dropna(subset=[x_var, metric])
            # Step 1: Fit a linear regression model
            X = working_data[[x_var]].values  # Independent variable (species richness)
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
            working_data['highlight'] = np.abs(working_data[f'{metric}_residuals']) > (2 * std_residual)

            # Print the cases with large differences
            # print("Cases with large differences (residuals):")
            # print(working_data[working_data['highlight']])

            working_data.to_csv(os.path.join(outpath, f'{metric}.csv'))
            plot_dist_of_metric(working_data, f'{metric}_residuals', out_path=os.path.join(outpath, 'dists', f'{metric}.jpg'))


def plot_regression_curves():
    # TODO: Save these as a grid
    # Optional: Plotting the results
    plt.scatter(working_data[x_var], working_data[metric], color='blue', label='Observed')
    plt.plot(working_data[x_var], working_data['expected_diversity'], color='red', label='Expected (Model Fit)')
    plt.title('Species Richness vs Diversity')
    plt.xlabel('Species Richness')
    plt.ylabel(metric)
    plt.legend()
    plt.savefig(os.path.join(outpath, f'{metric}.jpg'), dpi=300)
    plt.close()

    # Plot residuals
    plt.scatter(working_data[x_var], working_data[f'{metric}_residuals'], color='green', label='Residuals')
    plt.axhline(0, color='red', linestyle='--')
    plt.title('Residuals (Observed - Expected)')
    plt.xlabel('Species Richness')
    plt.ylabel('Residuals')
    plt.savefig(os.path.join(outpath, f'{metric}_residuals.jpg'), dpi=300)
    plt.close()


if __name__ == '__main__':
    # plot_distributions()
    # outliers()
    using_models()
