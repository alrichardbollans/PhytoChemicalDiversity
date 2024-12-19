import os

import pandas as pd
import statsmodels.api as sm
from pathlib import Path

from matplotlib import pyplot as plt
from scipy.stats import spearmanr
import seaborn as sns

from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import METRICS, RARE_METRICS



def get_working_data():
    phy_div_data = pd.read_csv(
        os.path.join('..', 'collect_and_compile_data', 'get_phylogenetic_diversities', 'outputs', 'group_data', 'native_regions_transformed.csv'),
        index_col=0)
    trait_data = pd.read_csv(
        os.path.join('..', 'collect_and_compile_data', 'get_diversity_metrics', 'outputs', 'group_data', 'native_regions_transformed.csv'),
        index_col=0)

    trait_data = trait_data.rename(columns={
        'Assigned_group': 'Group'
    })

    working_data = pd.merge(phy_div_data, trait_data, on='Group', how='inner')
    if 'number_of_species_in_data_and_tree' in working_data.columns:
        issues = working_data[working_data['number_of_species_in_group'] != working_data['number_of_species_in_data_and_tree']]

        if len(issues) > 0:
            print(issues)
            raise ValueError
    working_data = working_data.rename(columns={'phylogenetic_diversity':'Phylogenetic Diversity'})
    return working_data


# def get_group_data(metric: str, scale=True):
#     working_data = get_working_data()
#
#     # Define the independent variables (features) and dependent variable ('H')
#     X = working_data[used_vars]
#     y = working_data[metric]
#     if scale:
#
#         # Initialize the scaler
#         scaler = StandardScaler()
#
#         # Standardize the independent variables
#         X_scaled = scaler.fit_transform(X)
#
#         data = pd.DataFrame(X_scaled,
#                             index=X.index,
#                             columns=scaler.get_feature_names_out())
#     else:
#         data = X
#     data[metric] = y
#
#     Path(os.path.join('temp_outputs', metric)).mkdir(parents=True, exist_ok=True)
#     data.to_csv(os.path.join('temp_outputs', metric, 'native_regions.csv'))
#     return data


# def plot_relationships(data, metric, tag: str):
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#
#     # Create a figure with two subplots
#     plt.figure(figsize=(12, 6))
#
#     # Plot the relationship between PD and H
#
#     sns.regplot(x='Phylogenetic Diversity', y=metric, data=data, scatter_kws={'s': 20}, line_kws={'color': 'blue'})
#     # sns.regplot(x='Phylogenetic Diversity', y=metric, data=data, scatter=False, line_kws={'color': 'orange'}, order=3)
#
#     plt.ylabel(metric)
#     plt.xlabel('Phylogenetic Diversity')
#
#     # Show the plots
#     plt.tight_layout()
#     plt.savefig(os.path.join('outputs', 'reg_plots', f'{metric}_{tag}_phylogenetic_diversity_correlation.png'))
#
#     # Create a figure with two subplots
#     plt.figure(figsize=(12, 6))
#
#     # Plot the relationship between PD and H
#
#     sns.regplot(x='number_of_species_in_group', y=metric, data=data, scatter_kws={'s': 20}, line_kws={'color': 'blue'})
#
#     plt.ylabel(metric)
#     plt.xlabel('number_of_species_in_group')
#
#     # Show the plots
#     plt.tight_layout()
#     plt.savefig(os.path.join('outputs', 'reg_plots', f'{metric}_{tag}_number_of_species_in_group_correlation.png'))


def f_test(data, metric: str, tag: str, outpath:str):
    """
    The F-test of overall significance indicates whether your regression model provides a better fit than a model that contains no independent variables
    """
    import seaborn as sns
    data = data[['Phylogenetic Diversity', metric]].dropna(subset=['Phylogenetic Diversity', metric], how='any')
    reg_data = data[['Phylogenetic Diversity']]
    reg_data = sm.add_constant(reg_data) # Add a constant as when data is scaled, y!=0 when x=0
    y = data[metric]
    model = sm.OLS(y, reg_data).fit()
    out = model.summary()
    Path(os.path.join(outpath, 'ftests')).mkdir(parents=True, exist_ok=True)
    with open(os.path.join(outpath, 'ftests', f'ftest_{tag}_{metric}.csv'), 'w') as f:
        f.write(out.as_csv())
    f_value = model.fvalue
    p_value = model.f_pvalue

    df = pd.DataFrame([f_value, p_value], index=['f_value', 'p_value'], columns=[f'{tag}_{metric}'])
    return df


#
# def plot_var_reg_for_data(data, tag: str):
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#
#     # Create a figure with two subplots
#     plt.figure(figsize=(12, 6))
#
#     # Plot the relationship between PD and H
#
#     sns.regplot(x='PD', y='number_of_species_in_group', data=data, scatter_kws={'s': 20}, line_kws={'color': 'blue'})
#
#     # Plot the relationship between number_of_species_in_group and H
#
#     plt.ylabel('number_of_species_in_group')
#     plt.xlabel('PD')
#
#     # Show the plots
#     plt.tight_layout()
#     plt.savefig(os.path.join('outputs', f'{tag}_var_correlation.png'))


def partial_correlation_analysis(data, metric: str, tag: str, method='spearman'):
    '''
    Partial Correlation: This method is actually helpful in the presence of multicollinearity, as it measures the unique relationship between
    each predictor and metric while controlling for the other predictor.

    partial correlation measures the degree of association between two random variables, with the effect of a set of controlling random variables removed.

    https://pingouin-stats.org/build/html/generated/pingouin.partial_corr.html
    :param data: scaled data
    :param metric:
    :return:
    '''
    data = data.dropna(subset=['Phylogenetic Diversity', metric], how='any')
    # Correlation analysis
    # Univariate analyses showing correlations exist
    correlation_matrix = data[['Phylogenetic Diversity', 'number_of_species_in_group', metric]].corr(method=method)[metric]
    correlation_matrix = correlation_matrix.loc[['Phylogenetic Diversity', 'number_of_species_in_group']]
    correlation_matrix.columns = [f'{tag}_{metric}']
    print("\nCorrelation Matrix:")
    print(correlation_matrix)

    phydiv = spearmanr(data['Phylogenetic Diversity'], data[metric])
    taxdiv = spearmanr(data['number_of_species_in_group'], data[metric])

    # just check calculations are same as for plots
    last_computed_value = correlation_matrix.loc['Phylogenetic Diversity']
    assert round(phydiv.correlation, 10) == round(last_computed_value, 10)

    spearmanr_df = pd.DataFrame([phydiv.pvalue, taxdiv.pvalue], index=['Phylogenetic Diversity', 'number_of_species_in_group'],
                                columns=[f'{tag}_{metric}'])

    # Then attempt to disassociate one from the other
    import pingouin as pg
    div_importance = pg.partial_corr(data=data, x='Phylogenetic Diversity', y=metric, covar=['number_of_species_in_group'],
                                     method=method)
    div_importance.index = ['Phylogenetic Diversity']
    # num_importance = pg.partial_corr(data=data, x='number_of_species_in_group', y=metric, covar=['Phylogenetic Diversity'],
    #                                  method=method)
    # num_importance.index = ['number_of_species_in_group']
    pg_df = div_importance
    pg_df.columns = [f'{metric}' + c for c in pg_df.columns]

    pc_r_df = pg_df[[f'{metric}r']]
    pc_r_df.columns = [metric]
    pc_p_df = pg_df[[f'{metric}p-val']]

    return correlation_matrix, spearmanr_df, pg_df, pc_r_df, pc_p_df


def main(metrics, outpath):
    ftests = pd.DataFrame()
    correlations = pd.DataFrame()
    correlations_p = pd.DataFrame()
    native_regions_pg_dfs = pd.DataFrame()
    partial_correlation_r_dfs = pd.DataFrame()
    partial_correlation_p_dfs = pd.DataFrame()
    for metric in metrics:
        tag = 'native_regions'
        working_data = get_working_data()
        # unscaled_data = get_group_data(metric, scale=False)
        f_df = f_test(working_data, metric=metric, tag=tag, outpath=outpath)
        ftests = pd.concat([ftests, f_df], axis=1)

        # scaled_data = get_group_data(metric)
        correlation_matrix, spearmanr_df, pg_df, pc_r_df, pc_p_df = partial_correlation_analysis(working_data, metric=metric, tag=tag)
        correlations = pd.concat([correlations, correlation_matrix], axis=1)
        correlations_p = pd.concat([correlations_p, spearmanr_df], axis=1)
        native_regions_pg_dfs = pd.concat([native_regions_pg_dfs, pg_df], axis=1)
        partial_correlation_r_dfs = pd.concat([partial_correlation_r_dfs, pc_r_df], axis=1)
        partial_correlation_p_dfs = pd.concat([partial_correlation_p_dfs, pc_p_df], axis=1)
    ftests.to_csv(os.path.join(outpath, 'ftests', 'results.csv'))
    correlations.to_csv(os.path.join(outpath, 'correlations.csv'))

    fig, ax = plt.subplots(figsize=(5, 1.5))
    sns.heatmap(correlations.loc[['Phylogenetic Diversity']], cmap='inferno', annot=True, cbar=False)
    plt.yticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, 'correlation_heatmap.jpg'), dpi=300)
    plt.close()

    correlations_p.to_csv(os.path.join(outpath, 'correlations_p.csv'))

    native_regions_pg_dfs.to_csv(os.path.join(outpath, 'native_regions_partil_correlations.csv'))

    partial_correlation_p_dfs.to_csv(os.path.join(outpath, 'partil_correlations_p_values.csv'))
    partial_correlation_r_dfs.to_csv(os.path.join(outpath, 'partial_correlation_r_values.csv'))

    fig, ax = plt.subplots(figsize=(5, 1.5))
    sns.heatmap(partial_correlation_r_dfs, cmap='inferno', annot=True, cbar=False)
    plt.yticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(outpath, 'partial_correlation_heatmap.jpg'), dpi=300)
    plt.close()


if __name__ == '__main__':
    main(METRICS, os.path.join('outputs', 'correlations_with_pd'))
    main(RARE_METRICS, os.path.join('outputs', 'correlations_with_pd', 'rare'))
