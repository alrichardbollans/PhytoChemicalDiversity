import os

import pandas as pd


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
    working_data = working_data.rename(columns={'phylogenetic_diversity': 'PD', 'number_of_species_in_group': 'SR'})
    return working_data


def _multiple_loess(metric):
    ## not using this
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
