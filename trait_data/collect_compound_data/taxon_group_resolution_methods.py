import os
from typing import List, Any

import pandas as pd
from fontTools.subset import subset
from pandas import DataFrame
from phytochempy.data_compilation_utilities import get_pathway_version_resolved_at_taxon_level
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_accepted_columns

from trait_data.collect_compound_data import FAMILIES_OF_INTEREST, NP_PATHWAYS, COMPOUND_ID_COL, get_npclassifier_pathway_columns_in_df, \
    raw_all_taxa_compound_csv

_output_path = resource_filename(__name__, 'outputs')
all_genus_compound_csv = os.path.join(_output_path, 'all_genus_compound_data.csv')
all_species_compound_csv = os.path.join(_output_path, 'all_species_compound_data.csv')
genus_pathway_data_csv = os.path.join(_output_path, 'genus_level_pathway_data.csv')
genus_distinct_pathway_data_csv = os.path.join(_output_path, 'genus_level_distinct_pathway_data.csv')
species_in_study_csv = os.path.join(_output_path, 'species_in_study.csv')


def get_relevant_deduplicated_data(taxa_compound_data: pd.DataFrame, comp_id_col: str, taxon_grouping: str, families: List[str]) -> pd.DataFrame:
    """
    Removes records with unknown compound IDs or taxon name.

    Drop duplicate records (by compound ID and taxon name).

    Drop records not in study families.

    :param taxa_compound_data: A pandas DataFrame containing the taxa and compound data.
    :param comp_id_col: The name of the column in taxa_compound_data containing the compound IDs.
    :param taxon_grouping: The name of the column in taxa_compound_data containing the taxon grouping information.
    :param families: A list of strings representing the families to filter the data by.
    :return: A processed pandas DataFrame containing the filtered metabolite data.

    """
    from phytochempy.compound_properties import sanitize_filename

    # Remove records without necessary data, as well as duplicates
    taxa_compound_data = taxa_compound_data.dropna(subset=[comp_id_col, taxon_grouping], how='any')

    cleaned = taxa_compound_data.drop_duplicates(
        subset=[taxon_grouping, comp_id_col],
        keep='first')

    processed_metabolite_data = cleaned[cleaned[wcvp_accepted_columns['family']].isin(families)]
    processed_metabolite_data['NPclassif_pathway_results'] = processed_metabolite_data['NPclassif_pathway_results'].apply(sanitize_filename)

    pathway_cols = get_npclassifier_pathway_columns_in_df(processed_metabolite_data)
    pathways = []
    for p in pathway_cols:
        processed_metabolite_data[p] = processed_metabolite_data[p].apply(sanitize_filename)
        pathways += processed_metabolite_data[p].dropna().tolist()
    pathways = set(pathways)
    assert len(pathways) == 7
    print(pathways)
    return processed_metabolite_data


def separate_into_pathway(df: pd.DataFrame, pathway: str) -> tuple[DataFrame, Any]:
    """
    Separate the given DataFrame into two separate DataFrames based on a specified pathway.

    :param df: The input DataFrame.
    :param pathway: The pathway to separate the DataFrame by.
    :return: A tuple containing two separate DataFrames: positives and negatives.
    """
    df = df.dropna(subset=['NPclassif_pathway_results'])
    pathway_cols = get_npclassifier_pathway_columns_in_df(df)

    positives = pd.DataFrame()
    for col in pathway_cols:
        positives = pd.concat([positives, df[df[col] == pathway]])
    negatives = df[~df.index.isin(positives.index)]
    assert len(positives) + len(negatives) == len(df)
    problems = positives[~positives['NPclassif_pathway_results'].str.contains(pathway)]
    assert len(problems) == 0
    problems = negatives[negatives['NPclassif_pathway_results'].str.contains(pathway)]
    assert len(problems) == 0

    return positives, negatives


def add_pathway_information_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add pathway information columns to a DataFrame.

    :param df: The DataFrame to add pathway information columns to.
    :return: The DataFrame with pathway information columns added.
    """
    df = df.dropna(subset=['NPclassif_pathway_results'])
    original_length = len(df)

    for pathway in NP_PATHWAYS:
        relevant_paths, other_paths = separate_into_pathway(df, pathway)
        relevant_paths[pathway] = 1
        other_paths[pathway] = 0

        pathway_df = pd.concat([relevant_paths, other_paths])[[COMPOUND_ID_COL, pathway]]
        pathway_df = pathway_df.drop_duplicates(subset=[COMPOUND_ID_COL, pathway])
        amibiguous_duplicates = pathway_df[pathway_df[COMPOUND_ID_COL].duplicated(keep=False)]
        if len(amibiguous_duplicates) > 0:  # A check that comps with same ID have been assigned same class
            print(
                f'WARNING: Some ambiguity for pathway: {pathway}. This is likely due to differing smiles strings for same given {COMPOUND_ID_COL}.')
            print(amibiguous_duplicates)

            raise ValueError

        df = df.merge(pathway_df, how='left', on=COMPOUND_ID_COL)
    assert len(df) == original_length  # Check no additions from merge

    return df


def split_multiple_pathways_into_duplicate_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Resolve cases with multiple assinged compounds by separating into different rows

    :param df: A pandas DataFrame containing multiple pathways encoded as binary values.
    :return: A pandas DataFrame with duplicate rows created for each pathway that has multiple occurrences.
    """
    issues_to_resolve = df[df[NP_PATHWAYS].sum(axis=1) > 1]
    # List to store new rows
    new_rows = []

    # Iterate through each row in the dataframe
    for index, row in issues_to_resolve.iterrows():
        # Count the number of 1s in the row
        count_ones = row[NP_PATHWAYS].sum()

        # If more than one 1 is found
        if count_ones > 1:
            one_cols = []
            for col in NP_PATHWAYS:
                if row[col] == 1:
                    one_cols.append(col)
            # Iterate through each column in the row
            for col in one_cols:

                new_row = row.copy()
                for p in NP_PATHWAYS:
                    new_row[p] = 0
                new_row[col] = 1
                new_rows.append(new_row.T)

    # Add new rows to dataframe
    resolution_df = pd.DataFrame()._append(new_rows)
    assert len(resolution_df[resolution_df[NP_PATHWAYS].sum(axis=1) > 1]) == 0
    # Drop duplicate rows
    df = df[df[NP_PATHWAYS].sum(axis=1) < 2]

    out_df = pd.concat([df, resolution_df])

    return out_df


def get_genus_level_version_for_all_pathways(df: pd.DataFrame, taxon_grouping='Genus', use_distinct: bool = False) -> pd.DataFrame:
    ## Generate genus data for all pathways

    if use_distinct:
        new_df = split_multiple_pathways_into_duplicate_rows(df)
    else:
        new_df = df.copy()

    out_df = pd.DataFrame()
    out_df[taxon_grouping] = new_df[taxon_grouping].unique()
    original_length = len(out_df)

    for pathway in NP_PATHWAYS:
        genus_pathway_df = get_pathway_version_resolved_at_taxon_level(new_df, pathway, taxon_grouping_col=taxon_grouping)
        # genus_pathway_df = genus_pathway_df[
        #     ['Genus', 'identified_compounds_count', f'identified_{pathway}_count', f'mean_identified_as_{pathway}']]
        if 'identified_compounds_count' not in out_df.columns:
            out_df = pd.merge(out_df, genus_pathway_df, on=['Genus'], how='left')
        else:
            out_df = pd.merge(out_df, genus_pathway_df, on=['Genus', 'identified_compounds_count'], how='left')
    assert len(out_df) == original_length

    return out_df


def _species():
    my_df = pd.read_csv(raw_all_taxa_compound_csv, index_col=0)
    my_df = my_df.drop(columns=[wcvp_accepted_columns['name'],
                                wcvp_accepted_columns['name_w_author'],
                                wcvp_accepted_columns['rank'],
                                wcvp_accepted_columns['parent_name'],
                                wcvp_accepted_columns['species_ipni_id'],
                                ])  # Drop these as this is now a 'genus' dataset
    processed = get_relevant_deduplicated_data(my_df, COMPOUND_ID_COL, wcvp_accepted_columns['species'], FAMILIES_OF_INTEREST)

    processed_with_pathway_columns = add_pathway_information_columns(processed)

    processed_with_pathway_columns.to_csv(all_species_compound_csv)
    processed_with_pathway_columns.describe(include='all').to_csv(os.path.join(_output_path, 'all_species_compound_data_summary.csv'))

    species_in_study = processed_with_pathway_columns.drop_duplicates(subset=['accepted_species'], keep='first')[
        ['accepted_family', 'Genus', 'accepted_species', 'accepted_species_w_author']]
    issues = processed_with_pathway_columns[~processed_with_pathway_columns['accepted_species'].isin(processed['accepted_species'].values)]
    print(issues)
    species_in_study.to_csv(species_in_study_csv)


def resolve_compound_data_to_group(taxon_grouping: str):
    my_df = pd.read_csv(raw_all_taxa_compound_csv, index_col=0)
    my_df = my_df.drop(columns=[wcvp_accepted_columns['name'],
                                wcvp_accepted_columns['name_w_author'],
                                wcvp_accepted_columns['rank'],
                                wcvp_accepted_columns['species'],
                                wcvp_accepted_columns['species_w_author'],
                                wcvp_accepted_columns['parent_name'],
                                wcvp_accepted_columns['species_ipni_id'],
                                ])  # Drop these as this is now a 'genus' dataset
    processed = get_relevant_deduplicated_data(my_df, COMPOUND_ID_COL, taxon_grouping, FAMILIES_OF_INTEREST)

    group_compound_data = add_pathway_information_columns(processed)

    # Remove genera with only a single known compound prior to calculations
    counts = group_compound_data.value_counts(taxon_grouping)
    groups_with_single_compounds = pd.DataFrame({taxon_grouping: counts.index, 'N': counts.values})
    groups_with_single_compounds = groups_with_single_compounds[groups_with_single_compounds['N'] < 2][taxon_grouping].values.tolist()

    group_compound_data = group_compound_data[~group_compound_data['Genus'].isin(groups_with_single_compounds)]

    if taxon_grouping == 'Genus':
        group_compound_data.to_csv(all_genus_compound_csv)
        group_compound_data.describe(include='all').to_csv(os.path.join(_output_path, 'all_genus_compound_data_summary.csv'))
        ## Add a compound summary
        compound_summary = group_compound_data.drop_duplicates(subset=[COMPOUND_ID_COL])[[
            'example_compound_name', 'InChIKey', 'InChIKey_simp', 'Standard_SMILES', 'SMILES', 'CAS ID', 'NPclassif_class_results',
            'NPclassif_superclass_results',
            'NPclassif_pathway_results', 'NPclassif_isglycoside', 'Terpenoids', 'Fatty_acids', 'Polyketides', 'Carbohydrates',
            'Amino_acids_and_Peptides',
            'Shikimates_and_Phenylpropanoids', 'Alkaloids']]
        compound_summary.to_csv(os.path.join(_output_path, 'genus_compound_info.csv'))
        compound_summary[compound_summary['Polyketides'] == 1].reset_index(drop=True).to_csv(os.path.join(_output_path, 'Polyketides_info.csv'))
        compound_summary.describe(include='all').to_csv(os.path.join(_output_path, 'genus_compound_summary.csv'))

    if taxon_grouping == 'Genus':
        genus_pathway_data = get_genus_level_version_for_all_pathways(group_compound_data)
        genus_pathway_data.to_csv(genus_pathway_data_csv)
        genus_pathway_data.describe(include='all').to_csv(os.path.join(_output_path, 'genus_data_summary.csv'))
        import seaborn
        import matplotlib.pyplot as plt

        seaborn.displot(genus_pathway_data, x="identified_compounds_count", binwidth=1)
        plt.savefig(os.path.join('outputs', 'N_distribution.jpg'), dpi=300)
        plt.close()
    ## Get version where rows are repeated to account for single compounds with multiple pathways
    group_pathway_data = get_genus_level_version_for_all_pathways(group_compound_data, use_distinct=True)
    if taxon_grouping == 'Genus':
        group_pathway_data.to_csv(genus_distinct_pathway_data_csv)

    return group_compound_data, group_pathway_data


def _genera():
    resolve_compound_data_to_group('Genus')


if __name__ == '__main__':
    _genera()
    _species()
