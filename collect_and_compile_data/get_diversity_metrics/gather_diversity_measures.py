import os

import pandas as pd
from phytochempy.chemical_diversity_metrics import get_pathway_based_diversity_measures, calculate_FAD_measures, compile_rarified_calculations
from phytochempy.compound_properties import NP_PATHWAYS
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_accepted_columns

from collect_and_compile_data.collect_compound_data import all_species_compound_csv, \
    COMPOUND_ID_COL, FAMILIES_OF_INTEREST

FAD_INDICES = ['FAD', 'MFAD', 'APWD']
PATHWAY_INDICES = ['H', 'Hbc', 'G', 'J']
METRICS = PATHWAY_INDICES + FAD_INDICES
RARE_METRICS = [f'{x}_Rare' for x in METRICS]
_output_path = resource_filename(__name__, 'outputs')

max_number_of_pathways = len(NP_PATHWAYS)


def transform_compiled_data(compiled_data: pd.DataFrame, tag: str):
    from sklearn.preprocessing import PowerTransformer
    transformer = PowerTransformer(method='yeo-johnson')
    compiled_data = compiled_data.set_index('Assigned_group', drop=True)
    transformed_data = transformer.fit_transform(
        compiled_data[METRICS + RARE_METRICS + ['number_of_species_in_group', 'GroupSize_FAD', 'GroupSize_Pathways']])
    # Convert the transformed data back into a DataFrame
    df_transformed = pd.DataFrame(transformed_data,
                                  columns=METRICS + RARE_METRICS + ['number_of_species_in_group', 'GroupSize_FAD', 'GroupSize_Pathways'])
    df_transformed['Assigned_group'] = compiled_data.index
    df_transformed = df_transformed[
        ['Assigned_group', 'number_of_species_in_group'] + METRICS + RARE_METRICS + ['GroupSize_FAD', 'GroupSize_Pathways']]
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
    groups_with_enough_pathway_compounds = abundance_diversity[abundance_diversity['GroupSize_Pathways'] >= max_number_of_pathways][
        'Assigned_group'].unique().tolist()

    data_with_enough_compounds = working_data[working_data['Assigned_group'].isin(groups_with_enough_pathway_compounds)]
    FAD_measures = calculate_FAD_measures(data_with_enough_compounds, compound_grouping='Assigned_group')
    groups_with_enough_FAD_compounds = FAD_measures[FAD_measures['GroupSize_FAD'] >= max_number_of_pathways]['Assigned_group'].unique().tolist()

    study_groups = set(groups_with_enough_pathway_compounds).intersection(set(groups_with_enough_FAD_compounds))
    data_with_enough_compounds = working_data[working_data['Assigned_group'].isin(study_groups)]

    dropped_groups = working_data[~working_data['Assigned_group'].isin(study_groups)]['Assigned_group'].unique().tolist()
    original_groups = working_data['Assigned_group'].unique().tolist()
    assert set(dropped_groups).union(set(study_groups)) == set(original_groups)
    dropped_info = pd.DataFrame(
        [[str(dropped_groups), str(study_groups), str(original_groups)], [len(dropped_groups), len(study_groups), len(original_groups)]])
    dropped_info.to_csv(
        os.path.join('outputs', 'group_data', f'dropped_group_info_{tag}.csv'))

    fad_rare, pathway_rare = compile_rarified_calculations(data_with_enough_compounds, 'Assigned_group', max_number_of_pathways, COMPOUND_ID_COL,
                                                           1000)

    mean_values = mean_values[mean_values['Assigned_group'].isin(study_groups)]

    compiled_data = pd.merge(mean_values, abundance_diversity, how='left', on='Assigned_group', validate='one_to_one')
    compiled_data = pd.merge(compiled_data, FAD_measures, how='left', on='Assigned_group', validate='one_to_one')
    compiled_data = pd.merge(compiled_data, fad_rare, how='left', on='Assigned_group', validate='one_to_one')
    compiled_data = pd.merge(compiled_data, pathway_rare, how='left', on='Assigned_group', validate='one_to_one')

    for g in study_groups:
        assert g in compiled_data['Assigned_group'].values

    # compiled_data = pd.read_csv(os.path.join('outputs', 'group_data', f'{tag}.csv'), index_col=0)
    compiled_data = compiled_data.dropna(subset=METRICS, how='all')
    compiled_data.to_csv(os.path.join('outputs', 'group_data', f'{tag}.csv'))
    compiled_data.describe(include='all').to_csv(os.path.join('outputs', 'group_data', f'{tag}_summary.csv'))

    transform_compiled_data(compiled_data, tag)
