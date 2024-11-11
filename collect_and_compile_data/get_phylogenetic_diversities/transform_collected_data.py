import os.path

import pandas as pd


def transform_compiled_data(tag: str):
    compiled_data = pd.read_csv(os.path.join('outputs', 'group_data', f'{tag}.csv'))
    compiled_data.describe(include='all').to_csv(os.path.join('outputs', 'group_data', f'{tag}_summary.csv'))
    compiled_data = compiled_data.sort_values(by='Group')

    # First check the number of species match in the different datasets
    equivalent_trait_data = pd.read_csv(os.path.join('..', 'get_diversity_metrics', 'outputs', 'group_data', f'{tag}.csv'))
    compiled_data = compiled_data[compiled_data['Group'].isin(equivalent_trait_data['Assigned_group'].values)]
    equivalent_trait_data = equivalent_trait_data.set_index('Assigned_group')
    compiled_data = compiled_data.set_index('Group')
    pd.testing.assert_series_equal(equivalent_trait_data['number_of_species_in_group'], compiled_data['number_of_species_in_data_and_tree'],
                                   check_dtype=False, check_names=False)

    from sklearn.preprocessing import PowerTransformer
    transformer = PowerTransformer(method='yeo-johnson')

    transformed_data = transformer.fit_transform(compiled_data[["phylogenetic_diversity", "group_age"]])
    # Convert the transformed data back into a DataFrame
    df_transformed = pd.DataFrame(transformed_data, columns=["phylogenetic_diversity", "group_age"])
    df_transformed['Group'] = compiled_data.index
    df_transformed = df_transformed[['Group', "phylogenetic_diversity", "group_age"]]
    df_transformed.to_csv(os.path.join('outputs', 'group_data', f'{tag}_transformed.csv'))


if __name__ == '__main__':
    transform_compiled_data('native_regions')
