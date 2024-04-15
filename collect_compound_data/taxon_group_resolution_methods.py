import math
import os
from typing import List, Any

import pandas as pd
from pandas import DataFrame
from pkg_resources import resource_filename
from wcvp_download import wcvp_accepted_columns

from collect_compound_data import all_taxa_compound_csv, FAMILIES_OF_INTEREST, NP_PATHWAYS, COMPOUND_ID_COL, get_npclassifier_pathway_columns_in_df

_output_path = resource_filename(__name__, 'outputs')
genus_pathway_data_csv = os.path.join(_output_path, 'genus_level_pathway_data.csv')
genus_distinct_pathway_data_csv = os.path.join(_output_path, 'genus_level_distinct_pathway_data.csv')
processed_pathway_species_data_csv = os.path.join(_output_path, 'processed_with_pathway_columns.csv')


def get_relevant_deduplicated_data(taxa_compound_data: pd.DataFrame, comp_id_col: str, taxon_grouping: str, families: List[str]) -> pd.DataFrame:
    """
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


def get_genus_level_version_for_pathway(df: pd.DataFrame, pathway: str, taxon_grouping='Genus'):
    """
    :param df: A pandas DataFrame containing the data
    :param pathway: A string specifying the column name of the class
    :param taxon_grouping: A string specifying the taxon grouping column name (default is 'Genus')
    :return: A pandas DataFrame containing genus level version for pathway

    This method calculates the genus level version for a given pathway based on the provided DataFrame and parameters.
    It performs several transformations and computations on the DataFrame to derive the final result.
    """
    expected_mean = df[pathway].mean()

    genera_df = df.copy()

    num_genera_tested = len(genera_df[taxon_grouping].unique().tolist())

    # count labelled species
    counts = genera_df.value_counts(taxon_grouping)
    counts.name = 'identified_compounds_count'
    genera_df = pd.merge(genera_df, counts, how='left', left_on=taxon_grouping, right_index=True)
    N_class_col = f'identified_{pathway}_count'
    genera_df[N_class_col] = genera_df.groupby([taxon_grouping])[pathway].transform('sum')
    mean_col = f'mean_identified_as_{pathway}'
    genera_df[mean_col] = genera_df.groupby([taxon_grouping])[pathway].transform('mean')
    expected_mean_col = f'expected_total_mean_for_{pathway}'
    genera_df[expected_mean_col] = expected_mean

    # Normalised mean for some analysis following:
    # Daniele Micci-Barreca, ‘A Preprocessing Scheme for High-Cardinality Categorical Attributes in Classification and Prediction Problems’,
    # ACM SIGKDD Explorations Newsletter 3, no. 1 (July 2001): 27–32, https://doi.org/10.1145/507533.507538.
    # The means for each class are highly unreliable for small counts
    # We can use a blend of posterior and prior probabilties to improve this i.e. low evidence examples are corrected towards the population mean.
    # In below, k (as in original paper) determines half of the sample size for which we trust the mean estimate
    # f denotes the smoothing effect to balance categorical average vs prior. Higher value means stronger regularization.
    # Here we use the defaults used in the target encoder library
    def weighting_factor(given_val: float, k: int = 20, f: float = 10) -> float:
        denom = 1 + (math.e ** ((k - given_val) / f))
        return 1 / denom

    factor_col = f'weigthing_factor_for_{pathway}'
    genera_df[factor_col] = genera_df['identified_compounds_count'].apply(weighting_factor)

    norm_mean_col = f'norm_mean_identified_as_{pathway}'
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
        genus_pathway_df = get_genus_level_version_for_pathway(new_df, pathway, taxon_grouping=taxon_grouping)
        # genus_pathway_df = genus_pathway_df[
        #     ['Genus', 'identified_compounds_count', f'identified_{pathway}_count', f'mean_identified_as_{pathway}']]
        if 'identified_compounds_count' not in out_df.columns:
            out_df = pd.merge(out_df, genus_pathway_df, on=['Genus'], how='left')
        else:
            out_df = pd.merge(out_df, genus_pathway_df, on=['Genus', 'identified_compounds_count'], how='left')
    assert len(out_df) == original_length

    return out_df


if __name__ == '__main__':
    my_df = pd.read_csv(all_taxa_compound_csv, index_col=0)
    processed = get_relevant_deduplicated_data(my_df, 'SMILES', 'Genus', FAMILIES_OF_INTEREST)
    processed_with_pathway_columns = add_pathway_information_columns(processed)
    processed_with_pathway_columns.to_csv(processed_pathway_species_data_csv)
    processed_with_pathway_columns.describe(include='all').to_csv(os.path.join('outputs', 'processed_pathway_summary.csv'))
    genus_pathway_data = get_genus_level_version_for_all_pathways(processed_with_pathway_columns)
    genus_pathway_data.to_csv(genus_pathway_data_csv)
    genus_pathway_data.describe(include='all').to_csv(os.path.join(_output_path, 'genus_data_summary.csv'))
    import seaborn
    import matplotlib.pyplot as plt
    seaborn.displot(genus_pathway_data, x="identified_compounds_count", binwidth=1)
    plt.savefig(os.path.join('outputs', 'N_distribution.jpg'), dpi=300)
    plt.close()
    ## Get version where rows are repeated to account for single compounds with multiple pathways
    genus_pathway_data = get_genus_level_version_for_all_pathways(processed_with_pathway_columns, use_distinct=True)
    genus_pathway_data.to_csv(genus_distinct_pathway_data_csv)
