import os
from typing import List

import pandas as pd
from phytochempy.compound_properties import get_npclassifier_classes_from_df, \
    get_npclassifier_pathway_columns_in_df, NP_PATHWAYS
from phytochempy.data_compilation_utilities import merge_and_tidy_compound_datasets, tidy_final_dataset
from phytochempy.knapsack_searches import get_knapsack_data
from phytochempy.wikidata_searches import get_wikidata
from pkg_resources import resource_filename
from wcvpy.wcvp_download import wcvp_accepted_columns

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')
_tidied_outputs_folder = resource_filename(__name__, 'tidied_outputs')

_output_path = resource_filename(__name__, 'outputs')

raw_all_taxa_compound_csv = os.path.join(_temp_outputs_path, 'raw_all_taxa_compound_data.csv')

if not os.path.isdir(_temp_outputs_path):
    os.mkdir(_temp_outputs_path)
if not os.path.isdir(_tidied_outputs_folder):
    os.mkdir(_tidied_outputs_folder)
if not os.path.isdir(_output_path):
    os.mkdir(_output_path)

FAMILIES_OF_INTEREST = ['Gelsemiaceae', 'Gentianaceae', 'Apocynaceae', 'Loganiaceae', 'Rubiaceae']
WCVP_VERSION = '12'
COMPOUND_ID_COL = 'Standard_SMILES'
species_in_study_csv = os.path.join(_output_path, 'species_in_study.csv')
all_species_compound_csv = os.path.join(_output_path, 'all_species_compound_data.csv')


def collect_all_compound_data():
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
    tidy_final_dataset(with_npclass_classes, _tidied_outputs_folder, raw_all_taxa_compound_csv, COMPOUND_ID_COL)

    summary = pd.read_csv(raw_all_taxa_compound_csv, index_col=0)
    summary.describe(include='all').to_csv(os.path.join(_temp_outputs_path, 'raw_all_taxa_compound_data_summary.csv'))

    duplicate_smiles = summary[summary.duplicated(subset=['Standard_SMILES', wcvp_accepted_columns['name_w_author']], keep=False)]
    duplicate_smiles.sort_values(by='SMILES').to_csv('duplicate_smiles.csv')


def refine_to_species():
    my_df = pd.read_csv(raw_all_taxa_compound_csv, index_col=0)
    my_df = my_df.drop(columns=[wcvp_accepted_columns['name'],
                                wcvp_accepted_columns['name_w_author'],
                                wcvp_accepted_columns['rank'],
                                wcvp_accepted_columns['parent_name'],
                                wcvp_accepted_columns['species_ipni_id'],
                                ])  # Drop these as this is now a 'genus' dataset

    my_df  = my_df.dropna(subset=['NPclassif_pathway_results'])


    def get_relevant_deduplicated_data(taxa_compound_data: pd.DataFrame, comp_id_col: str, compound_grouping: str) -> pd.DataFrame:
        """
        Removes records with unknown compound IDs or taxon name.

        Drop duplicate records (by compound ID and taxon name).

        Drop records not in study families.

        :param taxa_compound_data: A pandas DataFrame containing the taxa and compound data.
        :param comp_id_col: The name of the column in taxa_compound_data containing the compound IDs.
        :param compound_grouping: The name of the column in taxa_compound_data containing the taxon grouping information.
        :param families: A list of strings representing the families to filter the data by.
        :return: A processed pandas DataFrame containing the filtered metabolite data.

        """
        from phytochempy.compound_properties import sanitize_filename

        # Remove records without necessary data, as well as duplicates
        taxa_compound_data = taxa_compound_data.dropna(subset=[comp_id_col, compound_grouping], how='any')

        cleaned = taxa_compound_data.drop_duplicates(
            subset=[compound_grouping, comp_id_col],
            keep='first')

        processed_metabolite_data = cleaned[cleaned[wcvp_accepted_columns['family']].isin(FAMILIES_OF_INTEREST)]
        # processed_metabolite_data['NPclassif_pathway_results'] = processed_metabolite_data['NPclassif_pathway_results'].apply(sanitize_filename)

        pathway_cols = get_npclassifier_pathway_columns_in_df(processed_metabolite_data)
        pathways = []
        for p in pathway_cols:
            # processed_metabolite_data[p] = processed_metabolite_data[p].apply(sanitize_filename)
            pathways += processed_metabolite_data[p].dropna().tolist()
        pathways = set(pathways)
        assert len(pathways) == 7
        print(pathways)
        return processed_metabolite_data

    processed = get_relevant_deduplicated_data(my_df, COMPOUND_ID_COL, wcvp_accepted_columns['species'])

    processed.to_csv(all_species_compound_csv)
    processed.describe(include='all').to_csv(os.path.join(_output_path, 'all_species_compound_data_summary.csv'))

    species_in_study = processed.drop_duplicates(subset=['accepted_species'], keep='first')[
        ['accepted_family', 'Genus', 'accepted_species', 'accepted_species_w_author']]
    issues = processed[~processed['accepted_species'].isin(species_in_study['accepted_species'].values)]
    if len(issues) > 0:
        print(issues)
    species_in_study.to_csv(species_in_study_csv)


if __name__ == '__main__':
    # collect_all_compound_data()
    refine_to_species()
