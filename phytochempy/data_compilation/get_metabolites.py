import os.path

import numpy as np
import pandas as pd
from wcvp_download import wcvp_accepted_columns
from wcvp_name_matching import get_genus_from_full_name

from phytochempy.compound_properties import CLASSYFIRE_OUTPUT_COLUMNS, get_classyfire_classes_from_df, get_compound_info_from_chembl_apm_assays, \
    simplify_inchi_key, add_chembl_data_to_compound_df, get_bioavailability_rules, add_bioavailability_rules_to_df
from phytochempy.knapsack_searches import get_knapsack_compounds_for_list_of_families
from phytochempy.wikidata_searches import generate_wikidata_search_query, submit_query, tidy_wikidata_output

NP_CLASSIFIER_COLUMNS = [
    'NPclassif_class_results', 'NPclassif_superclass_results',
    'NPclassif_pathway_results']
compound_class_columns = CLASSYFIRE_OUTPUT_COLUMNS + NP_CLASSIFIER_COLUMNS


def get_wikidata(wiki_data_id_for_order: str, temp_output_csv: str, tidied_output_csv: str, limit: int = 100000):
    # Example usage
    my_query = generate_wikidata_search_query(wiki_data_id_for_order, limit)
    submit_query(my_query, temp_output_csv)
    tidy_wikidata_output(temp_output_csv, tidied_output_csv)


def get_knapsack_data(families_of_interest: list, temp_output_path: str, tidied_output_csv: str):
    get_knapsack_compounds_for_list_of_families(families_of_interest, temp_output_path, tidied_output_csv)


def merge_and_tidy_compound_datasets(datasets: list, output_csv: str):
    all_metabolites_in_taxa = pd.concat(datasets)

    ### Format
    all_metabolites_in_taxa = all_metabolites_in_taxa.sort_values(
        by=wcvp_accepted_columns['name']).reset_index(drop=True)
    start_cols = ['accepted_name_w_author', 'Metabolite', 'InChIKey', 'SMILES', 'Molecular formula', 'CAS ID',
                  'chembl_id', 'Source']

    all_metabolites_in_taxa = all_metabolites_in_taxa[
        start_cols + [col for col in all_metabolites_in_taxa.columns if
                      col not in start_cols]]

    # Tidy final list
    chem_id_cols = ['SMILES', 'InChIKey', 'CAS ID']
    # Remove cases with no identifiers
    all_metabolites_in_taxa = all_metabolites_in_taxa.dropna(subset=chem_id_cols, how='all')

    ### Fill in matching details for compound ID columns from other dataframes
    def fill_match_ids(given_col: str):
        cols_to_do = [cl for cl in chem_id_cols if cl != given_col]
        # Step 1: Create a mapping of cols_to_do[0] and cols_to_do[1] to 'SMILES' values for non-empty rows
        mapping1 = {}
        mapping2 = {}
        for index, row in all_metabolites_in_taxa.iterrows():
            if row[given_col] == row[given_col]:
                if row[cols_to_do[0]] == row[cols_to_do[0]]:
                    mapping1[row[cols_to_do[0]]] = row[given_col]
                if row[cols_to_do[1]] == row[cols_to_do[1]]:
                    mapping2[row[cols_to_do[1]]] = row[given_col]

        if any(x in mapping1.keys() for x in ['', np.nan, 'nan']):
            raise ValueError
        if any(x in mapping2.keys() for x in ['', np.nan, 'nan']):
            raise ValueError

        # Step 2: Iterate through the DataFrame to fill in the empty 'SMILES' values
        for index, row in all_metabolites_in_taxa.iterrows():
            if row[given_col] != row[given_col]:
                m1 = row[cols_to_do[0]]
                m2 = row[cols_to_do[1]]
                # prioritise inchikey and smiles
                if m1 in mapping1:
                    v = mapping1[m1]
                    all_metabolites_in_taxa.at[index, given_col] = v
                elif m2 in mapping2:
                    v = mapping2[m2]
                    all_metabolites_in_taxa.at[index, given_col] = v

    for c_id in chem_id_cols:
        fill_match_ids(c_id)

    all_metabolites_in_taxa['InChIKey_simp'] = all_metabolites_in_taxa['InChIKey'].apply(simplify_inchi_key)
    all_metabolites_in_taxa.to_csv(output_csv)
    return all_metabolites_in_taxa


def add_chembl_data(df: pd.DataFrame, temp_csv: str, out_csv: str = None):
    get_compound_info_from_chembl_apm_assays(temp_csv)
    df_with_assay_data = add_chembl_data_to_compound_df(df, temp_csv, out_csv)
    return df_with_assay_data


def add_classyfire_info(df: pd.DataFrame, _temp_output_path: str, output_csv: str = None):
    ### Get classyfire info
    # Server was down as of 13/12/23
    all_metabolites_with_classyfire_info = get_classyfire_classes_from_df(df, 'SMILES',
                                                                          tempout_dir=_temp_output_path)

    if output_csv is not None:
        all_metabolites_with_classyfire_info.to_csv(output_csv)
    return all_metabolites_with_classyfire_info


def add_bioavailability_info(df: pd.DataFrame, out_csv: str = None):
    bio_av = add_bioavailability_rules_to_df(df, 'SMILES')

    if out_csv is not None:
        bio_av.to_csv(out_csv)
    return bio_av


def get_manual_files_to_upload(df: pd.DataFrame, _temp_output_path):
    ''' For data without, AFAIK, no API. This generates files to upload to obtain information'''

    ### Get NPclassifier info
    # https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22NPCLASSIFIER%22%7D
    # https://ccms-ucsd.github.io/GNPSDocumentation/api/#structure-natural-product-classifier-np-classifier
    all_metabolites_to_send_to_npclassifier = df[['SMILES']].dropna().drop_duplicates(
        keep='first')
    all_metabolites_to_send_to_npclassifier.to_csv(
        os.path.join(_temp_output_path, 'smiles_for_np_classifier.csv'))

    ### Get MAIP info
    # https://www.ebi.ac.uk/chembl/maip/
    all_metabolites_to_send_to_maip = df[['SMILES']].dropna().drop_duplicates(
        keep='first')
    all_metabolites_to_send_to_maip['id'] = all_metabolites_to_send_to_maip['SMILES']
    all_metabolites_to_send_to_maip.to_csv(os.path.join(_temp_output_path, 'smiles_for_MAIP.csv'))


def add_manual_info_files(df: pd.DataFrame, npclassifier_output_file: str, maip_output_file: str):
    ### Get NPclassifier info
    np_classif_results = pd.read_csv(npclassifier_output_file,
                                     sep='\t').drop_duplicates(keep='first')
    rename_dict = {'smiles': 'SMILES'}
    for c in np_classif_results.columns:
        if c != 'smiles':
            rename_dict[c] = 'NPclassif_' + c
    np_classif_results = np_classif_results.rename(columns=rename_dict).dropna(subset='SMILES')

    ### Get MAIP values
    maip_results = pd.read_csv(maip_output_file)
    maip_results = maip_results.rename(columns={'smiles': 'SMILES', 'model_score': 'MAIP_model_score'})
    maip_results = maip_results[['SMILES', 'MAIP_model_score']].dropna(subset='SMILES')

    all_metabolites_with_info = pd.merge(df, np_classif_results[~np_classif_results['SMILES'].isna()], how='left',
                                         on='SMILES')
    all_metabolites_with_info = pd.merge(all_metabolites_with_info, maip_results[~maip_results['SMILES'].isna()], how='left',
                                         on='SMILES')

    return all_metabolites_with_info


def tidy_and_check_final_dataset(pre_final_df: pd.DataFrame, _temp_output_path: str, final_taxa_compound_csv: str,
                                 final_compound_info_csv: str) -> None:
    ## Tidy data a bit

    # Add genus column
    pre_final_df['Genus'] = pre_final_df[wcvp_accepted_columns['name']].apply(
        get_genus_from_full_name)

    ### Checks
    def has_different_maipscore(group):
        return group['MAIP_model_score'].nunique() > 1

    # Check cases with the same simplified inchikey_simp value have same SMILEs/MAIP score. This is overwhelmingly the case. Similarly for compound classes.
    inchi_problems1 = pre_final_df.groupby('InChIKey_simp').filter(has_different_maipscore).drop_duplicates(
        subset=['InChIKey_simp', 'SMILES'],
        keep='first').sort_values(
        by='InChIKey_simp')

    if len(inchi_problems1.index) > 0:
        inchi_problems1.to_csv(os.path.join(_temp_output_path, 'same_inchisimp_diff_smiles.csv'))
        print(
            f'WARNING: Same some compounds with same InChIsimp and different smiles. See {os.path.join(_temp_output_path, "same_inchisimp_diff_smiles.csv")}')

    def drop_repeat_id(df, id_col):
        # Remove known duplicates according to ID
        df = df[df[id_col].isna() | (~df.duplicated(subset=[id_col, wcvp_accepted_columns['name_w_author']], keep='first'))]
        return df

    def drop_repeat_ids(df):
        # Drop cases without smiles or inchi values where CAS ID is repeated
        df = df[(~df['SMILES'].isna()) | (~df['InChIKey_simp'].isna()) | ~df[df['CAS ID'].notnull()].duplicated(
            subset=['CAS ID', wcvp_accepted_columns['name_w_author']], keep='first')]
        return df

    pre_final_df = pre_final_df.sort_values(by='SMILES').reset_index(drop=True)
    ## Remove InChIKey_simp repeats
    ## See APMCompoundScreen library for analysis of this
    pre_final_df = drop_repeat_id(pre_final_df, 'InChIKey')
    pre_final_df = drop_repeat_ids(pre_final_df)
    taxa_with_all_info = pre_final_df.sort_values(by=[wcvp_accepted_columns['name_w_author'], 'Metabolite']).reset_index(drop=True)

    taxa_with_all_info.to_csv(final_taxa_compound_csv)
    # Note that sometimes classifications are the same but strings given in different order.
    output_compound_info = ['InChIKey_simp', 'active_chembl_compound'] + compound_class_columns + ['MAIP_model_score', 'lipinski_pass', 'veber_pass']

    # Remove repeated inchisimps and keep the most active/those with activity info
    taxa_with_all_info[
        output_compound_info].sort_values(by=['active_chembl_compound'], ascending=False).drop_duplicates(subset=['InChIKey_simp'],
                                                                                                          keep='first').reset_index(drop=True).to_csv(
        final_compound_info_csv)


if __name__ == '__main__':
    ### Example workflow

    # Define context
    families = ['Apocynaceae', 'Loganiaceae', 'Rubiaceae', ]
    wiki_data_id_for_order = 'Q21754'
    temp_outputs_folder = 'temp'
    tidied_outputs_folder = 'tidied'

    ## Get compound-taxa pair data
    # get_wikidata(wiki_data_id_for_order, os.path.join(temp_outputs_folder, 'wikidata.csv'), os.path.join(tidied_outputs_folder, 'wikidata.csv'))
    # get_knapsack_data(families, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'knapsack_data.csv'))

    ## Merge and tidy the data
    tidy_wiki_data = pd.read_csv(os.path.join(tidied_outputs_folder, 'wikidata.csv'))
    tidy_knapsack_data = pd.read_csv(os.path.join(tidied_outputs_folder, 'knapsack_data.csv'))
    all_compounds_in_taxa = merge_and_tidy_compound_datasets([tidy_wiki_data, tidy_knapsack_data],
                                                             os.path.join(tidied_outputs_folder, 'merged_data.csv'))

    ## Add extra information related to the compound properties
    # These steps can be included/removed as needed
    # For the longer processes, to avoid repeats you can simply read the associated temp_output if the step has already been run
    with_chembl_data = add_chembl_data(all_compounds_in_taxa, os.path.join(temp_outputs_folder, 'chembl.csv'))
    # with_chembl_data = pd.read_csv(os.path.join(temp_outputs_folder, 'chembl.csv'), index_col=0)
    with_classyfire_classes = add_classyfire_info(with_chembl_data, temp_outputs_folder, os.path.join(tidied_outputs_folder, 'classyfire.csv'))
    with_bioavailibility = add_bioavailability_info(with_classyfire_classes, os.path.join(tidied_outputs_folder, 'bioavailibility.csv'))

    get_manual_files_to_upload(with_bioavailibility, temp_outputs_folder)
    all_info = add_manual_info_files(with_bioavailibility, 'example_np_classifier_file.tsv', 'example_maip_file.csv')

    ### Then tidy and output final dataset
    tidy_and_check_final_dataset(all_info, tidied_outputs_folder, os.path.join('outputs', 'all_taxa_compound_data.csv'),
                                 os.path.join('outputs', 'all_compound_data.csv'))
