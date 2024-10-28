import os

import pandas as pd
from phytochempy.chemical_diversity_metrics import get_pathway_based_diversity_measures, calculate_FAD_measures
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_accepted_columns

from collect_and_compile_data.collect_compound_data import all_species_compound_csv, \
    COMPOUND_ID_COL, FAMILIES_OF_INTEREST

FAD_INDICES = ['FAD', 'MFAD', 'APWD']
PATHWAY_INDICES = ['H', 'Hbc', 'G']
METRICS = PATHWAY_INDICES + FAD_INDICES
_output_path = resource_filename(__name__, 'outputs')


def transform_compiled_data(compiled_data: pd.DataFrame, tag: str):
    from sklearn.preprocessing import PowerTransformer
    transformer = PowerTransformer(method='yeo-johnson')
    compiled_data = compiled_data.set_index('Assigned_group', drop=True)
    transformed_data = transformer.fit_transform(compiled_data[METRICS + ['number_of_species_in_group', 'N']])
    # Convert the transformed data back into a DataFrame
    df_transformed = pd.DataFrame(transformed_data, columns=METRICS + ['number_of_species_in_group', 'N'])
    df_transformed['Assigned_group'] = compiled_data.index
    df_transformed = df_transformed[['Assigned_group', 'number_of_species_in_group'] + METRICS + ['N']]
    df_transformed.to_csv(os.path.join('outputs', 'group_data', f'{tag}_transformed.csv'))


def resolve_traits_to_group(df: pd.DataFrame, tag: str):
    # resolve traits to group
    assert len(df[df.duplicated(subset=['accepted_species', 'Assigned_group'])].index) == 0

    df['number_of_species_in_group'] = df[['Assigned_group', 'accepted_species']].groupby('Assigned_group').transform('count')
    df.to_csv(os.path.join('outputs', 'group_info', f'{tag}.csv'))

    mean_values = df[['Assigned_group', 'number_of_species_in_group']].groupby(
        'Assigned_group').mean()
    mean_values = mean_values.reset_index()

    def check_means(x):
        if x != int(x):
            raise ValueError
        else:
            pass

    mean_values['number_of_species_in_group'].apply(check_means)
    print(mean_values)

    # After mean values have been calculated, add compound data
    compound_data = pd.read_csv(all_species_compound_csv, index_col=0)
    working_data = pd.merge(compound_data, df, how='left', on='accepted_species', validate='many_to_many')
    working_data = working_data[working_data[wcvp_accepted_columns['family']].isin(FAMILIES_OF_INTEREST)]
    working_data = working_data.dropna(subset='Assigned_group')

    abundance_diversity = get_pathway_based_diversity_measures(working_data, 'Assigned_group', COMPOUND_ID_COL)
    FAD_measures = calculate_FAD_measures(working_data, taxon_grouping='Assigned_group')

    compiled_data = pd.merge(mean_values, abundance_diversity, how='left', on='Assigned_group', validate='one_to_one')
    compiled_data = pd.merge(compiled_data, FAD_measures, how='left', on='Assigned_group', validate='one_to_one')

    for g in list(set(working_data['Assigned_group'].values.tolist())) + list(set(FAD_measures['Assigned_group'].values.tolist())):
        assert g in compiled_data['Assigned_group'].values

    compiled_data.to_csv(os.path.join('outputs', 'group_data', f'{tag}.csv'))
    # compiled_data = pd.read_csv(os.path.join('outputs', 'group_data', f'{tag}.csv'))
    transform_compiled_data(compiled_data, tag)


def genera():
    my_df = pd.read_csv(all_species_compound_csv, index_col=0)[['accepted_species', 'Genus']].drop_duplicates()
    my_df['Assigned_group'] = my_df['Genus']
    resolve_traits_to_group(my_df, tag='Genus')


if __name__ == '__main__':
    genera()
