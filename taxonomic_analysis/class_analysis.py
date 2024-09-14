import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from pkg_resources import resource_filename
from statsmodels.stats.proportion import proportion_confint

from trait_data.collect_compound_data import NP_PATHWAYS, genus_pathway_data_csv

_output_path = resource_filename(__name__, 'outputs')


def binom_mean_ci(n, p, alpha):
    """
    Calculate the confidence interval for the mean of a binomial distribution using clopper pearson.

    Parameters:
    n (int): Number of trials
    p (float): Probability of success
    alpha (float): Significance level

    Returns:
    tuple: Lower and upper bounds of the confidence interval
    """
    n_success = n * p
    ci = proportion_confint(n_success, n, alpha, method='beta')

    return ci


def get_CI_df(genera_df: pd.DataFrame, pathway: str):
    genera_pathway_df = genera_df[
        ['Genus', 'identified_compounds_count', f'identified_{pathway}_count', f'mean_identified_as_{pathway}', f'expected_total_mean_for_{pathway}']]

    population_mean = genera_pathway_df[f'expected_total_mean_for_{pathway}'].tolist()[0]
    max_number_tested = genera_pathway_df['identified_compounds_count'].max()

    ninety_nine_lower_bounds = {}
    ninety_nine_upper_bounds = {}
    for i in range(1, max_number_tested + 1):
        low, upp = binom_mean_ci(i, population_mean, alpha=0.01)

        ninety_nine_lower_bounds[i] = low

        ninety_nine_upper_bounds[i] = upp
    lower_series = pd.Series(list(ninety_nine_lower_bounds.values()), index=ninety_nine_lower_bounds.keys())
    lower_series.name = 'lower_bound'
    CI_df = pd.DataFrame(lower_series)
    CI_df['count'] = ninety_nine_lower_bounds.keys()
    CI_df['upper_bound'] = pd.Series(list(ninety_nine_upper_bounds.values()), index=ninety_nine_upper_bounds.keys())
    return CI_df, genera_pathway_df


def plot_CI_df(ci_df: pd.DataFrame, genera_pathway_df: pd.DataFrame, pathway: str, log=False):
    '''
    Funnel plot as in milliken_plants_2021.
    Also note approximations of distribution and confidence intervals are less useful for small values of n and means near extremes i.e. if np(1 − p)<10.
    '''
    population_mean = genera_pathway_df[f'expected_total_mean_for_{pathway}'].tolist()[0]

    if log:
        line_var = 'logarithm_count'
        ci_df[line_var] = np.log10(ci_df['count'])
    else:
        line_var = 'count'

    sns.lineplot(x=line_var, y='lower_bound', data=ci_df, linestyle='--', label='Lower bound')
    splot = sns.lineplot(x=line_var, y='upper_bound', data=ci_df, linestyle='--', label='Upper bound')
    splot.axhline(population_mean, linestyle='--', label='Population Mean', color='black')
    if log:
        scatter_var = 'logarithm_identified_compounds_count'

        genera_pathway_df['logarithm_identified_compounds_count'] = np.log10(genera_pathway_df['identified_compounds_count'])
    else:
        scatter_var = 'identified_compounds_count'
    merged = pd.merge(ci_df, genera_pathway_df, left_on='count', right_on='identified_compounds_count', how='right')
    overactives = merged[merged[f'mean_identified_as_{pathway}'] > merged['upper_bound']]  # ['Genus'].tolist()
    underactives = merged[merged[f'mean_identified_as_{pathway}'] < merged['lower_bound']]  # ['Genus'].tolist()

    overactive_genera = overactives['Genus'].tolist()
    underactive_genera = underactives['Genus'].tolist()
    print(f'Mean: {pathway}')
    print('Overactives: ', overactive_genera)
    print('underactives: ', underactive_genera)
    # Add scatter plot
    sns.scatterplot(x=scatter_var, y=f'mean_identified_as_{pathway}', data=genera_pathway_df, color='grey')
    overaactive_colour = '#d95f02'
    sns.scatterplot(x=scatter_var, y=f'mean_identified_as_{pathway}', data=overactives, color=overaactive_colour)

    sns.scatterplot(x=scatter_var, y=f'mean_identified_as_{pathway}', data=underactives, color='#7570b3')

    # Annotate specific points
    for i, row in genera_pathway_df.iterrows():
        genus = row['Genus']
        if genus in overactive_genera or genus in underactive_genera:  # Replace with the genera you want to annotate
            x = genera_pathway_df.loc[i, scatter_var]
            y = genera_pathway_df.loc[i, f'mean_identified_as_{pathway}']
            if log:
                threshold = 2
            else:
                threshold = 100
            if x > threshold:
                if x < 0:
                    # Calculate the annotation offset based on the values
                    xytext_x = np.random.uniform(-50, 50)  # Adjust the range as needed
                    xytext_y = np.random.uniform(-50, 50)  # Adjust the range as needed

                    arrow_props = dict(facecolor='black')

                    plt.annotate(genus, (x, y), xytext=(xytext_x, xytext_y), textcoords='offset points', arrowprops=arrow_props)
                else:
                    plt.annotate(genus, (x, y), xytext=(5, 5), textcoords='offset points')
    if log:
        plt.xlabel('Log10 Number of Identified Compounds')
    else:
        plt.xlabel('Number of Identified Compounds')

    plt.ylabel(f'Ratio of {pathway.capitalize()}s')
    plt.ylim(-0.1, 1.1)
    # plt.title('Mean Species Activity')
    plt.legend()
    if log:
        plt.savefig(os.path.join(_output_path, f'{pathway}_confidence_intervals_log.jpg'), dpi=300)
    else:
        plt.savefig(os.path.join(_output_path, f'{pathway}_confidence_intervals.jpg'), dpi=300)

    plt.close()

def notable_genera():
    given_genera_df = pd.read_csv(genus_pathway_data_csv, index_col=0)

    # Top 5 for each index
    for pathway in NP_PATHWAYS:
        top_5 = given_genera_df.sort_values(by=f'mean_identified_as_{pathway}', ascending=False).head(5)
        top_5.to_csv(os.path.join(_output_path, 'top_5_' + pathway + '.csv'))

        top_5 = given_genera_df.sort_values(by=f'norm_mean_identified_as_{pathway}', ascending=False).head(5)
        top_5.to_csv(os.path.join(_output_path, 'top_5_' + pathway + '_normed.csv'))

def main():
    # Note with current set up this uses data where genera with single compounds (N=1) have been removed.

    notable_genera()
    given_genera_df = pd.read_csv(genus_pathway_data_csv, index_col=0)
    for pway in NP_PATHWAYS:
        my_ci_df, my_genera_pathway_df = get_CI_df(given_genera_df, pway)
        plot_CI_df(my_ci_df, my_genera_pathway_df, pway)
        plot_CI_df(my_ci_df, my_genera_pathway_df, pway, log=True)


if __name__ == '__main__':
    main()
