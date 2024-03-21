import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pkg_resources import resource_filename

from library_info_and_data_import import CLASSES_OF_INTEREST, processed_compounds_output_path

_output_path = resource_filename(__name__, 'outputs')

from scipy.stats import norm


def binom_mean_ci(n, p, alpha):
    """
    Calculate the confidence interval for the mean of a binomial distribution using the normal approximation.

    Parameters:
    n (int): Number of trials
    p (float): Probability of success
    alpha (float): Significance level

    Returns:
    tuple: Lower and upper bounds of the confidence interval
    """

    std_dev = np.sqrt((p * (1 - p) / n))

    z_crit = norm.ppf(1 - alpha / 2)
    margin_of_error = z_crit * std_dev

    lower_bound = p - margin_of_error
    upper_bound = p + margin_of_error

    return lower_bound, upper_bound


def get_CI_df(genera_df, expected_mean):
    max_number_tested = genera_df['tested_compounds_count'].max()

    ninety_nine_lower_bounds = {}
    ninety_nine_upper_bounds = {}
    for i in range(1, max_number_tested + 1):
        low, upp = binom_mean_ci(i, expected_mean, alpha=0.01)

        ninety_nine_lower_bounds[i] = low

        ninety_nine_upper_bounds[i] = upp
    lower_series = pd.Series(list(ninety_nine_lower_bounds.values()), index=ninety_nine_lower_bounds.keys())
    lower_series.name = 'lower_bound'
    CI_df = pd.DataFrame(lower_series)
    CI_df['count'] = ninety_nine_lower_bounds.keys()
    CI_df['upper_bound'] = pd.Series(list(ninety_nine_upper_bounds.values()), index=ninety_nine_upper_bounds.keys())
    return CI_df, genera_df


def plot_CI_df(ci_df, genera_df, comp_class):
    '''
    Funnel plot as in milliken_plants_2021.
    Also note approximations of distribution and confidence intervals are less useful for small values of n e.g. <5.
    '''

    sns.lineplot(x='count', y='lower_bound', data=ci_df, linestyle='--')
    sns.lineplot(x='count', y='upper_bound', data=ci_df, linestyle='--')

    merged = pd.merge(ci_df, genera_df, left_on='count', right_on='tested_compounds_count', how='right')
    overactives = merged[merged['mean_compound_value'] > merged['upper_bound']]['Genus'].tolist()
    underactives = merged[merged['mean_compound_value'] < merged['lower_bound']]['Genus'].tolist()
    print('Overactives: ', overactives)
    print('underactives: ', underactives)
    # Add scatter plot
    sns.scatterplot(x='tested_compounds_count', y='mean_compound_value', data=genera_df, label=f'Ratio of {comp_class.capitalize()}s', color='#d95f02')

    # Annotate specific points
    for i, row in genera_df.iterrows():
        genus = row['Genus']
        if genus in overactives or genus in underactives:  # Replace with the genera you want to annotate
            x = genera_df.loc[i, 'tested_compounds_count']
            y = genera_df.loc[i, 'mean_compound_value']
            if x < 0:
                # Calculate the annotation offset based on the values
                xytext_x = np.random.uniform(-50, 50)  # Adjust the range as needed
                xytext_y = np.random.uniform(-50, 50)  # Adjust the range as needed

                arrow_props = dict(facecolor='black')

                plt.annotate(genus, (x, y), xytext=(xytext_x, xytext_y), textcoords='offset points', arrowprops=arrow_props)
            else:
                plt.annotate(genus, (x, y), xytext=(5, 5), textcoords='offset points')

    plt.xlabel('Number of Known Compounds')
    plt.ylabel(f'Ratio of Known Compounds Which Are {comp_class.capitalize()}s')
    plt.ylim(-0.1, 1.1)
    # plt.title('Mean Species Activity')
    plt.legend()
    plt.savefig(os.path.join(_output_path, f'{comp_class}_confidence_intervals.jpg'), dpi=300)
    plt.close()


def main():
    for comp_class in CLASSES_OF_INTEREST:
        genera_df = pd.read_csv(os.path.join(processed_compounds_output_path, f'genus_level_{comp_class}_information.csv'))
        activity_mean = genera_df['expected_total_mean'].tolist()[0]
        ci_df, genera_df = get_CI_df(genera_df, activity_mean)
        plot_CI_df(ci_df, genera_df,comp_class)


if __name__ == '__main__':
    main()
