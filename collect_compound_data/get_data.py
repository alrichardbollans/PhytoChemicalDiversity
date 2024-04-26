import os
from typing import List

import pandas as pd
from phytochempy.compound_properties import get_npclassifier_result_columns_in_df, get_npclassifier_classes_from_df
from phytochempy.data_compilation_utilities import merge_and_tidy_compound_datasets, tidy_final_dataset
from phytochempy.knapsack_searches import get_knapsack_data
from phytochempy.wikidata_searches import get_wikidata
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_accepted_columns

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')
_tidied_outputs_folder = resource_filename(__name__, 'tidied_outputs')

_output_path = resource_filename(__name__, 'outputs')

all_taxa_compound_csv = os.path.join(_output_path, 'all_taxa_compound_data.csv')

if not os.path.isdir(_temp_outputs_path):
    os.mkdir(_temp_outputs_path)
if not os.path.isdir(_tidied_outputs_folder):
    os.mkdir(_tidied_outputs_folder)
if not os.path.isdir(_output_path):
    os.mkdir(_output_path)

FAMILIES_OF_INTEREST = ['Gelsemiaceae', 'Gentianaceae', 'Apocynaceae', 'Loganiaceae', 'Rubiaceae']
COMPOUND_ID_COL = 'Standard_SMILES'
NP_PATHWAYS = ['Terpenoids', 'Fatty_acids', 'Polyketides', 'Carbohydrates', 'Amino_acids_and_Peptides', 'Shikimates_and_Phenylpropanoids',
               'Alkaloids']


def get_npclassifier_pathway_columns_in_df(df: pd.DataFrame) -> List[str]:
    """
    :param df: A pandas DataFrame containing NPclassifier result columns.
    :return: A list of pathway columns in the given DataFrame.
    """
    pathway_cols = get_npclassifier_result_columns_in_df(df)
    cols = []
    for p in pathway_cols:
        if 'pathway' in p and p != 'NPclassif_pathway_results':
            cols.append(p)
    return cols


if __name__ == '__main__':
    # Define context
    wiki_data_id_for_order = 'Q21754'

    ## Get compound-taxa pair data
    # get_wikidata(wiki_data_id_for_order, os.path.join(_temp_outputs_path, 'wikidata.csv'), os.path.join(_tidied_outputs_folder, 'wikidata.csv'))
    # get_knapsack_data(FAMILIES_OF_INTEREST, _temp_outputs_path, os.path.join(_tidied_outputs_folder, 'knapsack_data.csv'))

    ## Merge and tidy the data
    tidy_wiki_data = pd.read_csv(os.path.join(_tidied_outputs_folder, 'wikidata.csv'), index_col=0)
    tidy_knapsack_data = pd.read_csv(os.path.join(_tidied_outputs_folder, 'knapsack_data.csv'), index_col=0)
    all_compounds_in_taxa = merge_and_tidy_compound_datasets([tidy_wiki_data, tidy_knapsack_data],
                                                             os.path.join(_tidied_outputs_folder, 'merged_data.csv'))

    ## Add extra information related to the compound properties

    # These steps can be included/removed as needed
    # For the longer processes, to avoid repeats you can simply read the associated temp_output if the step has already been run
    # all_compounds_in_taxa = pd.read_csv(os.path.join(_tidied_outputs_folder, 'merged_data.csv'), index_col=0)
    with_npclass_classes = get_npclassifier_classes_from_df(all_compounds_in_taxa, 'Standard_SMILES', _temp_outputs_path)
    with_npclass_classes.to_csv(os.path.join(_tidied_outputs_folder, 'npclassifier.csv'))
    pway_cols = get_npclassifier_pathway_columns_in_df(with_npclass_classes)

    ### Then tidy and output final dataset
    tidy_final_dataset(with_npclass_classes, _tidied_outputs_folder, all_taxa_compound_csv, COMPOUND_ID_COL)

    summary = pd.read_csv(all_taxa_compound_csv, index_col=0)
    summary.describe(include='all').to_csv(os.path.join(_output_path, 'all_taxa_compound_data_summary.csv'))

    duplicate_smiles = summary[summary.duplicated(subset=['Standard_SMILES', wcvp_accepted_columns['name_w_author']], keep=False)]
    duplicate_smiles.sort_values(by='SMILES').to_csv('duplicate_smiles.csv')