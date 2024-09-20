import os

import pandas as pd
import statsmodels.api as sm
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler

used_vars = ['phylogenetic_diversity', 'number_of_species_in_data_and_tree', 'Bio1', 'Animal Richness']


def get_group_data(tag: str, metric: str = 'H', scale=True):
    phy_div_data = pd.read_csv(os.path.join('..', 'trait_data', 'get_phylogeny', 'species_phylogeny', 'outputs', 'group_data', f'{str(tag)}.csv'),
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
    data.to_csv(os.path.join('temp_outputs', metric, f'{tag}.csv'))
    return data


def get_genera_data(metric: str = 'H', scale=True):
    phy_div_data = pd.read_csv(os.path.join('..', 'trait_data', 'get_phylogeny', 'species_phylogeny', 'outputs', 'phylogenetic_diversities.csv'),
                               index_col=0)  # only 158 becuase of number of singletons
    trait_data = pd.read_csv(os.path.join('..', 'trait_data', 'compile_trait_data', 'outputs', 'genus_trait_data.csv'), index_col=0)[
        ['Genus', 'num_species_in_data', metric, 'Bio1', 'Animal Richness']]

    working_data = pd.merge(phy_div_data, trait_data, on='Genus', how='inner')
    issues = working_data[working_data['num_species_in_data'] != working_data['number_of_species_in_data_and_tree']]

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

    from pathlib import Path
    Path(os.path.join('temp_outputs', metric)).mkdir(exist_ok=True)
    data.to_csv(os.path.join('temp_outputs', metric, f'genus_group.csv'))
    return data


def plot_relationships(data, metric):
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Create a figure with two subplots
    plt.figure(figsize=(12, 6))

    # Plot the relationship between phylogenetic_diversity and H

    sns.regplot(x='phylogenetic_diversity', y=metric, data=data, scatter_kws={'s': 20}, line_kws={'color': 'blue'})

    # Plot the relationship between number_of_species_in_data_and_tree and H

    sns.regplot(x='number_of_species_in_data_and_tree', y=metric, data=data, scatter_kws={'s': 20}, line_kws={'color': 'orange'})

    plt.ylabel(metric)
    plt.xlabel('')
    plt.legend()

    # Show the plots
    plt.tight_layout()
    plt.show()


#
# def dominance_analysis(it: int, metric: str = 'H'):
#     '''
#     Relative Importance Analysis: The dominance analysis we perform is designed to handle multicollinearity better than standard regression
#     coefficients, which is why it's valuable in this case.
#     :param it:
#     :param metric:
#     :return:
#     '''
#     ## See https://github.com/dominance-analysis/dominance-analysis/issues/33 for fix
#     X_scaled, y, data = get_group_data(it)
#     from dominance_analysis import Dominance
#     import plotly
#     dominance_regression = Dominance(data=data, target=metric, objective=1)
#     incr_variable_rsquare = dominance_regression.incremental_rsquare()
#     fig = dominance_regression.plot_incremental_rsquare()
#     plt.show()
#     out = dominance_regression.dominance_stats()
#     level = dominance_regression.dominance_level()
#     print(out)
#     print(level)

def f_test(data, metric: str, tag: str):
    """
    The F-test of overall significance indicates whether your regression model provides a better fit than a model that contains no independent variables
    """

    import seaborn as sns
    reg_data = data[['phylogenetic_diversity', 'number_of_species_in_data_and_tree']]
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


def partial_correlation_analysis(data, metric: str, tag: str, method='spearman'):
    '''
    Partial Correlation: This method is actually helpful in the presence of multicollinearity, as it measures the unique relationship between
    each predictor and H while controlling for the other predictor.

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
    num_importance = pg.partial_corr(data=data, x='number_of_species_in_data_and_tree', y=metric, covar=['phylogenetic_diversity'],
                                     method=method)
    num_importance.index = ['number_of_species_in_data_and_tree']
    pg_df = pd.concat([div_importance, num_importance], axis=0)
    pg_df.columns = [f'{tag}_{metric}'+c for c in pg_df.columns]


    # This just important for 'regions'
    # div_importance = pg.partial_corr(data=data, x='phylogenetic_diversity', y=metric, covar=['Bio1', 'Animal Richness'],
    #                        method=method)
    # num_importance = pg.partial_corr(data=data, x='number_of_species_in_data_and_tree', y=metric, covar=['Bio1', 'Animal Richness'],
    #                        method=method)
    # print(div_importance)
    # print(num_importance)

    # plot_relationships(data, metric)

    return correlation_matrix, spearmanr_df, pg_df


def main():

    metrics = ['FAD', 'MFAD', 'APWD', 'H', 'Hbc', 'G', 'J']
    ftests = pd.DataFrame()
    correlations = pd.DataFrame()
    correlations_p = pd.DataFrame()
    pg_dfs = pd.DataFrame()
    for metric in metrics:
        unscaled_genera_data = get_genera_data(scale=False, metric=metric)
        f_df = f_test(unscaled_genera_data, metric=metric, tag='Genera')
        ftests = pd.concat([ftests, f_df], axis=1)

        scaled_genera_data = get_genera_data(metric=metric)
        correlation_matrix, spearmanr_df, pg_df = partial_correlation_analysis(scaled_genera_data, metric=metric, tag='Genera')
        correlations = pd.concat([correlations, correlation_matrix], axis=1)
        correlations_p = pd.concat([correlations_p, spearmanr_df], axis=1)
        pg_dfs = pd.concat([pg_dfs, pg_df], axis=1)

        for tag in ['random_genera', 'random_regions', 'native_regions']:
            unscaled_data = get_group_data(tag, scale=False, metric=metric)
            f_df = f_test(unscaled_data, metric=metric, tag=tag)
            ftests = pd.concat([ftests, f_df], axis=1)
            scaled_data = get_group_data(tag, metric=metric)
            correlation_matrix, spearmanr_df, pg_df = partial_correlation_analysis(scaled_data, metric=metric, tag=tag)
            correlations = pd.concat([correlations, correlation_matrix], axis=1)
            correlations_p = pd.concat([correlations_p, spearmanr_df], axis=1)
            pg_dfs = pd.concat([pg_dfs, pg_df], axis=1)


    ftests.to_csv(os.path.join('outputs', 'ftests', 'results.csv'))
    correlations.to_csv(os.path.join('outputs', 'correlations.csv'))
    correlations_p.to_csv(os.path.join('outputs', 'correlations_p.csv'))
    pg_dfs.to_csv(os.path.join('outputs', 'partil_correlations.csv'))


if __name__ == '__main__':
    main()
