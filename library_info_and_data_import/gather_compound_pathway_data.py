import math
import os
import re

import numpy as np
import pandas as pd
from metabolite_properties import filter_rows_containing_compound_keyword
from metabolites import all_taxa_metabolites_csv, compound_class_columns
from pkg_resources import resource_filename
from wcvp_download import wcvp_accepted_columns

_inputs_path = resource_filename(__name__, 'inputs')

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')

processed_compounds_output_path = resource_filename(__name__, 'outputs')
_chembl_assay_summary_csv = os.path.join(processed_compounds_output_path, 'chembl_assay_summary.csv')
_processed_metabolite_csv = os.path.join(processed_compounds_output_path, 'processed_metabolites.csv')
_class_output_path = os.path.join(processed_compounds_output_path, 'classes')
_pathway_output_path = os.path.join(processed_compounds_output_path, 'pathways')
processed_metabolite_with_classes_csv = os.path.join(_class_output_path, 'processed_metabolites_with_classes.csv')
processed_metabolite_with_pathways_csv = os.path.join(_pathway_output_path, 'processed_metabolites_with_pathways.csv')
# Uniqueness is up to 'SMILES'
# TODO: Check usage of full key when collecting data
COMPOUND_ID_COL = 'SMILES'
FAMILIES_OF_INTEREST = ['Apocynaceae', 'Loganiaceae', 'Rubiaceae']
CLASSES_OF_INTEREST = ['alkaloid', 'steroid', 'terpenoid', 'peptide', 'lipid', 'polyketide']
MINOR_CLASSES = ['pyrrolizidine', 'quinoline']
output_compound_info = [wcvp_accepted_columns['species_w_author'], 'Genus',
                        'Metabolite', 'InChIKey', 'SMILES',
                        ] + compound_class_columns + [wcvp_accepted_columns['species'], wcvp_accepted_columns['family']]


def sanitize_filename(pathway: str, replace_spaces: bool = True):
    """
    Sanitize given pathway names as they are currently not very nice for handling filenames and/or importing exporting in R etc..

    Args:
        pathway (str): The filename to sanitize.
        replace_spaces (bool): Whether to replace spaces with underscores.

    Returns:
        str: The sanitized filename.
    """
    if pathway == pathway:
        # Remove leading and trailing whitespace
        pathway = pathway.strip()

        # Replace spaces with underscores (optional)
        if replace_spaces:
            pathway = pathway.replace(" ", "_")

        # Remove illegal characters
        sanitized_filename = re.sub(r'[\\/*?:"<>|]', "", pathway)

        # Limit filename length (optional, for example, to 255 characters)
        sanitized_filename = sanitized_filename[:255]

        return sanitized_filename
    else:
        return pathway


def get_processed_metabolite_data(comp_id_col: str, taxon_grouping: str = 'Genus'):
    # Get plant compound data
    all_taxa_metabolite_data = pd.read_csv(all_taxa_metabolites_csv, index_col=0)
    # And do some cleaning
    all_taxa_metabolite_data = all_taxa_metabolite_data.dropna(subset=[comp_id_col, taxon_grouping], how='any')

    cleaned = all_taxa_metabolite_data.drop_duplicates(
        subset=[taxon_grouping, comp_id_col],
        keep='first')
    processed_metabolite_data = cleaned[output_compound_info]
    processed_metabolite_data = processed_metabolite_data[processed_metabolite_data[wcvp_accepted_columns['family']].isin(FAMILIES_OF_INTEREST)]

    processed_metabolite_data['NPclassif_pathway_results'] = processed_metabolite_data['NPclassif_pathway_results'].apply(sanitize_filename)
    processed_metabolite_data['classyfire_superclass'] = processed_metabolite_data['classyfire_superclass'].apply(sanitize_filename)

    processed_metabolite_data.to_csv(_processed_metabolite_csv)

    class_suggestions = processed_metabolite_data['NPclassif_pathway_results'].unique()
    print('NPclassif_pathway_results class suggestions')
    print(class_suggestions)

    class_suggestions = processed_metabolite_data['classyfire_superclass'].unique()
    print('classyfire_superclass class suggestions')
    print(class_suggestions)
    return processed_metabolite_data


def add_class_information_columns(metabolite_data: pd.DataFrame):
    ''' Allows finding specific classes of interest

    :param metabolite_data:
    :return:
    '''
    original_length = len(metabolite_data)
    for comp_class in CLASSES_OF_INTEREST + MINOR_CLASSES:
        alks, non_alks, nans = filter_rows_containing_compound_keyword(metabolite_data, compound_class_columns,
                                                                       comp_class)
        assert len(nans) == 0
        alks[comp_class] = 1
        non_alks[comp_class] = 0

        class_df = pd.concat([alks, non_alks])[[COMPOUND_ID_COL, comp_class]]
        class_df = class_df.drop_duplicates(subset=[COMPOUND_ID_COL, comp_class])
        amibiguous_duplicates = class_df[class_df[COMPOUND_ID_COL].duplicated(keep=False)]
        if not len(amibiguous_duplicates) == 0:  # A check that comps with same ID have been assigned same class
            print(f'WARNING: Some ambiguity for: {comp_class}. This is likely due to differing smiles strings for same given Inchikey.')
            print('Positive case is given preference.')
            print(amibiguous_duplicates)
            if len(amibiguous_duplicates) > 50:
                raise ValueError
            class_df = class_df.drop_duplicates(subset=[COMPOUND_ID_COL])
        metabolite_data = metabolite_data.merge(class_df, how='left', on=COMPOUND_ID_COL)
        assert len(metabolite_data) == original_length  # Check no additions from merge
    metabolite_data.to_csv(processed_metabolite_with_classes_csv)
    metabolite_data.describe(include='all').to_csv(os.path.join(_class_output_path, 'processed_metabolite_class_summary.csv'))

    return metabolite_data


def add_pathway_information_columns(metabolite_data: pd.DataFrame, pathway_column: str = 'classyfire_superclass'):
    ''' More general method for pathways in data'''
    original_length = len(metabolite_data)
    nan_pathways = metabolite_data[metabolite_data[pathway_column].isna()]
    non_nan_pathways = metabolite_data[~metabolite_data[pathway_column].isna()]
    for pathway in metabolite_data[pathway_column].unique():
        relevant_paths = non_nan_pathways[non_nan_pathways[pathway_column] == pathway]
        other_paths = non_nan_pathways[non_nan_pathways[pathway_column] != pathway]
        relevant_paths[pathway] = 1
        other_paths[pathway] = 0
        nan_pathways[pathway] = np.nan

        class_df = pd.concat([relevant_paths, other_paths, nan_pathways])[[COMPOUND_ID_COL, pathway]]
        class_df = class_df.drop_duplicates(subset=[COMPOUND_ID_COL, pathway])
        amibiguous_duplicates = class_df[class_df[COMPOUND_ID_COL].duplicated(keep=False)]
        if not len(amibiguous_duplicates) == 0:  # A check that comps with same ID have been assigned same class
            print(f'WARNING: Some ambiguity for pathway: {pathway}. This is likely due to differing smiles strings for same given Inchikey.')
            print('Positive case is given preference.')
            print(amibiguous_duplicates)
            if len(amibiguous_duplicates) > 50:
                raise ValueError
            class_df = class_df.drop_duplicates(subset=[COMPOUND_ID_COL])
        metabolite_data = metabolite_data.merge(class_df, how='left', on=COMPOUND_ID_COL)
        assert len(metabolite_data) == original_length  # Check no additions from merge
    metabolite_data.to_csv(processed_metabolite_with_pathways_csv)
    metabolite_data.describe(include='all').to_csv(os.path.join(processed_compounds_output_path, 'processed_metabolite_pathway_summary.csv'))


def get_genus_level_version_for_compound_class(metabolite_data: pd.DataFrame, comp_class: str, taxon_grouping='Genus'):
    expected_mean = metabolite_data[comp_class].mean()

    genera_df = metabolite_data.copy()

    num_genera_tested = len(genera_df[taxon_grouping].unique().tolist())

    # count labelled species
    counts = genera_df.value_counts(taxon_grouping)
    counts.name = 'identified_compounds_count'
    genera_df = pd.merge(genera_df, counts, how='left', left_on=taxon_grouping, right_index=True)
    N_class_col = f'identified_{comp_class}_count'
    genera_df[N_class_col] = genera_df.groupby([taxon_grouping])[comp_class].transform('sum')
    mean_col = f'mean_identified_as_{comp_class}'
    genera_df[mean_col] = genera_df.groupby([taxon_grouping])[comp_class].transform('mean')
    expected_mean_col = f'expected_total_mean_for_{comp_class}'
    genera_df[expected_mean_col] = expected_mean

    # Normalised mean following:
    # Daniele Micci-Barreca, ‘A Preprocessing Scheme for High-Cardinality Categorical Attributes in Classification and Prediction Problems’,
    # ACM SIGKDD Explorations Newsletter 3, no. 1 (July 2001): 27–32, https://doi.org/10.1145/507533.507538.
    # The means for each class are highly unreliable for small counts
    # We can use a blend of posterior and prior probabilties to improve this.
    # In below, k (as in original paper) determines half of the sample size for which we trust the mean estimate
    # f denotes the smoothing effect to balance categorical average vs prior. Higher value means stronger regularization.
    # Here we use the defaults used in the target encoder library
    def weighting_factor(given_val: float, k: int = 20, f: float = 10) -> float:
        denom = 1 + (math.e ** ((k - given_val) / f))
        return 1 / denom

    factor_col = f'weigthing_factor_for_{comp_class}'
    genera_df[factor_col] = genera_df['identified_compounds_count'].apply(weighting_factor)

    norm_mean_col = f'norm_mean_identified_as_{comp_class}'
    genera_df[norm_mean_col] = (genera_df[factor_col] * genera_df[mean_col]) + (
            1 - genera_df[factor_col]) * expected_mean

    genera_df = genera_df[
        [taxon_grouping, 'identified_compounds_count', N_class_col, mean_col, expected_mean_col, factor_col,
         norm_mean_col]]
    genera_df = genera_df.reset_index(drop=True)

    genera_df = genera_df.drop_duplicates(subset=[taxon_grouping])
    genera_df.reset_index(drop=True, inplace=True)
    assert len(genera_df) == num_genera_tested
    return genera_df


def get_genus_level_version_for_all_classes(metabolite_data: pd.DataFrame, pathways: list, out_dir: str, taxon_grouping='Genus'):
    ## Generate genus data for all pathways
    out_df = pd.DataFrame()
    out_df[taxon_grouping] = metabolite_data[taxon_grouping].unique()
    original_length = len(out_df)
    for pathway in pathways:
        genus_pathway_df = get_genus_level_version_for_compound_class(metabolite_data, pathway, taxon_grouping=taxon_grouping)
        genus_pathway_df = genus_pathway_df[
            ['Genus', 'identified_compounds_count', f'identified_{pathway}_count', f'mean_identified_as_{pathway}']]
        if 'identified_compounds_count' not in out_df.columns:
            out_df = pd.merge(out_df, genus_pathway_df, on=['Genus'], how='left')
        else:
            out_df = pd.merge(out_df, genus_pathway_df, on=['Genus', 'identified_compounds_count'], how='left')
    assert len(out_df) == original_length

    out_df.to_csv(os.path.join(out_dir, 'genus_level_pathway_data.csv'))
    return out_df


def get_pathway_categories(pathway_column='classyfire_superclass'):
    metabolite_data = pd.read_csv(_processed_metabolite_csv, index_col=0)
    return metabolite_data[pathway_column].unique().tolist()


def main():
    get_processed_metabolite_data(COMPOUND_ID_COL)
    # add_class_information_columns(metabolite_data=pd.read_csv(_processed_metabolite_csv, index_col=0))
    #
    # df_with_class_columns = pd.read_csv(processed_metabolite_with_classes_csv, index_col=0)
    # for comp_class in CLASSES_OF_INTEREST + MINOR_CLASSES:
    #     get_genus_level_version_for_compound_class(df_with_class_columns, comp_class, _class_output_path)

    metabolite_data = pd.read_csv(_processed_metabolite_csv, index_col=0)
    metabolite_data = metabolite_data.dropna(subset=['classyfire_superclass'])
    add_pathway_information_columns(metabolite_data, pathway_column='classyfire_superclass')
    df_with_pathway_columns = pd.read_csv(processed_metabolite_with_pathways_csv, index_col=0)
    get_genus_level_version_for_all_classes(df_with_pathway_columns, df_with_pathway_columns['classyfire_superclass'].unique(), _pathway_output_path)


if __name__ == '__main__':
    main()
