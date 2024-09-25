import os

import pandas as pd
import statsmodels.api as sm
from pathlib import Path
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler

used_vars = ['phylogenetic_diversity', 'number_of_species_in_data_and_tree']


def get_group_data(tag: str, metric: str = 'H', scale=True):
    phy_div_data = pd.read_csv(os.path.join('..', 'trait_data', 'get_phylogenetic_diversities', 'outputs', 'group_data', f'{str(tag)}.csv'),
                               index_col=0)
    trait_data = pd.read_csv(os.path.join('..', 'trait_data', 'other_group_traits', 'outputs', 'group_data', f'{str(tag)}.csv'), index_col=0)

    trait_data = trait_data.rename(columns={
        'Assigned_group': 'Group'
    })

    working_data = pd.merge(phy_div_data, trait_data, on='Group', how='inner')
    issues = working_data[working_data['number_of_species_in_group'] != working_data['number_of_species_in_data_and_tree']]

    if len(issues) > 0:
        print(issues)
        raise ValueError

    # Define the independent variables (features) and dependent variable ('H')
    X = working_data[used_vars]
    y = working_data[metric]
    if scale:

        # Initialize the scaler
        scaler = StandardScaler()

        # Standardize the independent variables
        X_scaled = scaler.fit_transform(X)

        data = pd.DataFrame(X_scaled,
                            index=X.index,
                            columns=scaler.get_feature_names_out())
    else:
        data = X
    data[metric] = y

    Path(os.path.join('temp_outputs', metric)).mkdir(parents=True, exist_ok=True)
    data.to_csv(os.path.join('temp_outputs', metric, f'{tag}.csv'))
    return data


def plot_relationships(data, metric, tag: str):
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Create a figure with two subplots
    plt.figure(figsize=(12, 6))

    # Plot the relationship between phylogenetic_diversity and H

    sns.regplot(x='phylogenetic_diversity', y=metric, data=data, scatter_kws={'s': 20}, line_kws={'color': 'blue'})
    # sns.regplot(x='phylogenetic_diversity', y=metric, data=data, scatter=False, line_kws={'color': 'orange'}, order=3)

    plt.ylabel(metric)
    plt.xlabel('phylogenetic_diversity')

    # Show the plots
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', 'reg_plots', f'{metric}_{tag}_phylogenetic_diversity_correlation.png'))

    # Create a figure with two subplots
    plt.figure(figsize=(12, 6))

    # Plot the relationship between phylogenetic_diversity and H

    sns.regplot(x='number_of_species_in_data_and_tree', y=metric, data=data, scatter_kws={'s': 20}, line_kws={'color': 'blue'})

    plt.ylabel(metric)
    plt.xlabel('number_of_species_in_data_and_tree')

    # Show the plots
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', 'reg_plots', f'{metric}_{tag}_number_of_species_in_data_and_tree_correlation.png'))


def f_test(data, metric: str, tag: str):
    """
    The F-test of overall significance indicates whether your regression model provides a better fit than a model that contains no independent variables
    """

    import seaborn as sns
    reg_data = data[['phylogenetic_diversity']]
    # reg_data = sm.add_constant(reg_data) # Don't add a constant as we know that y=0 when X1,X2=0 when data isnt scaled (this is the case for all metrics)
    y = data[metric]
    model = sm.OLS(y, reg_data).fit()
    out = model.summary()
    with open(os.path.join('outputs', 'ftests', f'ftest_{tag}_{metric}.csv'), 'w') as f:
        f.write(out.as_csv())
    f_value = model.fvalue
    p_value = model.f_pvalue

    df = pd.DataFrame([f_value, p_value], index=['f_value', 'p_value'], columns=[f'{tag}_{metric}'])
    return df


def plot_var_reg_for_data(data, tag: str):
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Create a figure with two subplots
    plt.figure(figsize=(12, 6))

    # Plot the relationship between phylogenetic_diversity and H

    sns.regplot(x='phylogenetic_diversity', y='number_of_species_in_data_and_tree', data=data, scatter_kws={'s': 20}, line_kws={'color': 'blue'})

    # Plot the relationship between number_of_species_in_data_and_tree and H

    plt.ylabel('number_of_species_in_data_and_tree')
    plt.xlabel('phylogenetic_diversity')

    # Show the plots
    plt.tight_layout()
    plt.savefig(os.path.join('outputs', f'{tag}_var_correlation.png'))


def partial_correlation_analysis(data, metric: str, tag: str, method='spearman'):
    '''
    Partial Correlation: This method is actually helpful in the presence of multicollinearity, as it measures the unique relationship between
    each predictor and H while controlling for the other predictor.

    partial correlation measures the degree of association between two random variables, with the effect of a set of controlling random variables removed.

    https://pingouin-stats.org/build/html/generated/pingouin.partial_corr.html
    :param data: scaled data
    :param metric:
    :return:
    '''

    # Correlation analysis
    # Univariate analyses showing correlations exist
    correlation_matrix = data.corr(method=method)[[metric]]
    correlation_matrix = correlation_matrix.loc[['phylogenetic_diversity', 'number_of_species_in_data_and_tree']]
    correlation_matrix.columns = [f'{tag}_{metric}']
    print("\nCorrelation Matrix:")
    print(correlation_matrix)

    phydiv = spearmanr(data['phylogenetic_diversity'], data[metric])
    taxdiv = spearmanr(data['number_of_species_in_data_and_tree'], data[metric])
    spearmanr_df = pd.DataFrame([phydiv.pvalue, taxdiv.pvalue], index=['phylogenetic_diversity', 'number_of_species_in_data_and_tree'],
                                columns=[f'{tag}_{metric}'])

    # Then attempt to disassociate one from the other
    import pingouin as pg
    div_importance = pg.partial_corr(data=data, x='phylogenetic_diversity', y=metric, covar=['number_of_species_in_data_and_tree'],
                                     method=method)
    div_importance.index = ['phylogenetic_diversity']
    # num_importance = pg.partial_corr(data=data, x='number_of_species_in_data_and_tree', y=metric, covar=['phylogenetic_diversity'],
    #                                  method=method)
    # num_importance.index = ['number_of_species_in_data_and_tree']
    pg_df = div_importance
    pg_df.columns = [f'{tag}_{metric}' + c for c in pg_df.columns]


    return correlation_matrix, spearmanr_df, pg_df


def main():
    metrics = ['FAD', 'MFAD', 'APWD', 'H', 'Hbc', 'G', 'J']
    ftests = pd.DataFrame()
    correlations = pd.DataFrame()
    correlations_p = pd.DataFrame()
    native_regions_pg_dfs = pd.DataFrame()
    random_regions_pg_dfs = pd.DataFrame()
    for metric in metrics:

        for tag in ['native_regions', 'random_regions']:
            unscaled_data = get_group_data(tag, scale=False, metric=metric)
            f_df = f_test(unscaled_data, metric=metric, tag=tag)
            ftests = pd.concat([ftests, f_df], axis=1)
            plot_var_reg_for_data(unscaled_data, tag)
            plot_relationships(unscaled_data, metric, tag)

            scaled_data = get_group_data(tag, metric=metric)
            correlation_matrix, spearmanr_df, pg_df = partial_correlation_analysis(scaled_data, metric=metric, tag=tag)
            correlations = pd.concat([correlations, correlation_matrix], axis=1)
            correlations_p = pd.concat([correlations_p, spearmanr_df], axis=1)
            if tag == 'native_regions':
                native_regions_pg_dfs = pd.concat([native_regions_pg_dfs, pg_df], axis=1)
            elif tag == 'random_regions':
                random_regions_pg_dfs = pd.concat([random_regions_pg_dfs, pg_df], axis=1)
            else:
                raise ValueError
    ftests.to_csv(os.path.join('outputs', 'ftests', 'results.csv'))
    correlations.to_csv(os.path.join('outputs', 'correlations.csv'))
    correlations_p.to_csv(os.path.join('outputs', 'correlations_p.csv'))

    def sort_pg_data(pg_results):
        num_species_more_important = []
        phyl_div_more_important = []

        significant = []
        not_significant = []

    random_regions_pg_dfs.to_csv(os.path.join('outputs', 'random_partil_correlations.csv'))
    native_regions_pg_dfs.to_csv(os.path.join('outputs', 'native_regions_partil_correlations.csv'))


if __name__ == '__main__':
    main()
