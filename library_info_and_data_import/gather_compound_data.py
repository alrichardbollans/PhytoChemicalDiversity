import os

import pandas as pd
from metabolite_properties import filter_rows_containing_compound_keyword
from metabolites import all_taxa_metabolites_csv, compound_class_columns, chembl_apm_assay_info_csv
from pkg_resources import resource_filename
from wcvp_download import wcvp_accepted_columns

_inputs_path = resource_filename(__name__, 'inputs')

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')

processed_compounds_output_path = resource_filename(__name__, 'outputs')
_chembl_assay_summary_csv = os.path.join(processed_compounds_output_path, 'chembl_assay_summary.csv')
_processed_metabolite_csv = os.path.join(processed_compounds_output_path, 'processed_metabolites.csv')
processed_metabolite_with_classes_csv = os.path.join(processed_compounds_output_path, 'processed_metabolites_with_classes.csv')
# Uniqueness is up to 'InChIKey_simp' i.e. the connectivity information
# 'active_chembl_compound' values are preserved within COMPOUND_ID_COL -- has been checked and is built in this way
# TODO:Try with full key
COMPOUND_ID_COL = 'InChIKey'
FAMILIES_OF_INTEREST = ['Apocynaceae', 'Loganiaceae', 'Rubiaceae']
CLASSES_OF_INTEREST = ['alkaloid', 'steroid', 'quinine', 'terpenoid']

output_compound_info = [wcvp_accepted_columns['species_w_author'], 'Genus',
                        'example_metabolite_name', 'InChIKey', 'SMILES', COMPOUND_ID_COL, 'active_chembl_compound',
                        'mean_ic50_μM', 'lipinski_pass', 'veber_pass',
                        'MAIP_model_score',
                        ] + compound_class_columns + [wcvp_accepted_columns['species'], wcvp_accepted_columns['family']]


def get_processed_metabolite_data():
    ## First get assay data
    chembl_compound_assays = pd.read_csv(os.path.join(chembl_apm_assay_info_csv),
                                         index_col=0).dropna(subset=['InChIKey', 'assay_pchembl_value'], how='any')

    chembl_compound_assays_summary = chembl_compound_assays[['InChIKey_simp', 'mean_ic50_μM']].drop_duplicates(keep='first')
    assert len(chembl_compound_assays_summary['InChIKey_simp'].unique().tolist()) == len(chembl_compound_assays['InChIKey_simp'].unique().tolist())

    chembl_compound_assays_summary = chembl_compound_assays_summary.sort_values(by='InChIKey_simp').reset_index(drop=True)
    chembl_compound_assays_summary.to_csv(_chembl_assay_summary_csv)

    # Then get plant compound data
    all_taxa_metabolite_data = pd.read_csv(all_taxa_metabolites_csv, index_col=0)
    # And do some cleaning
    all_taxa_metabolite_data = all_taxa_metabolite_data.dropna(subset=[COMPOUND_ID_COL, wcvp_accepted_columns['species_w_author']], how='any')
    all_taxa_metabolite_data = all_taxa_metabolite_data.rename(columns={'Metabolite': 'example_metabolite_name'})

    cleaned = all_taxa_metabolite_data.drop_duplicates(
        subset=[wcvp_accepted_columns['species_w_author'], COMPOUND_ID_COL],
        keep='first')

    # Add mean bioassay results
    processed_metabolite_data = pd.merge(cleaned, chembl_compound_assays_summary, how='left',
                                         on=COMPOUND_ID_COL)
    processed_metabolite_data = processed_metabolite_data[output_compound_info]
    processed_metabolite_data = processed_metabolite_data[processed_metabolite_data[wcvp_accepted_columns['family']].isin(FAMILIES_OF_INTEREST)]
    processed_metabolite_data.to_csv(_processed_metabolite_csv)


def add_class_information():
    metabolite_data = pd.read_csv(_processed_metabolite_csv, index_col=0)
    original_length = len(metabolite_data)
    for comp_class in CLASSES_OF_INTEREST:
        alks, non_alks, nans = filter_rows_containing_compound_keyword(metabolite_data, compound_class_columns,
                                                                       comp_class)
        alks[comp_class] = 1
        non_alks[comp_class] = 0

        class_df = pd.concat([alks, non_alks])[[COMPOUND_ID_COL, comp_class]]
        class_df = class_df.drop_duplicates(subset=[COMPOUND_ID_COL, comp_class])
        if not len(class_df[class_df[COMPOUND_ID_COL].duplicated()]) == 0:  # A check that comps with same ID have been assigned same class
            print(comp_class)
            print(class_df[class_df[COMPOUND_ID_COL].duplicated(keep=False)])
        metabolite_data = metabolite_data.merge(class_df, how='left', on=COMPOUND_ID_COL)
        assert len(metabolite_data) == original_length  # Check no additions from merge
    metabolite_data.to_csv(processed_metabolite_with_classes_csv)
    metabolite_data.describe(include='all').to_csv(os.path.join(processed_compounds_output_path, 'processed_metabolite_summary.csv'))


def get_genus_level_version(comp_class: str):
    metabolite_data = pd.read_csv(processed_metabolite_with_classes_csv, index_col=0)
    expected_mean = metabolite_data[comp_class].mean()

    genera_df = metabolite_data.copy()

    num_genera_tested = len(genera_df['Genus'].unique().tolist())

    # count labelled species
    counts = genera_df.value_counts('Genus')
    counts.name = 'tested_compounds_count'
    genera_df = pd.merge(genera_df, counts, how='left', left_on='Genus', right_index=True)
    genera_df['hit_compounds_count'] = genera_df.groupby(['Genus'])[comp_class].transform('sum')
    genera_df['mean_compound_value'] = genera_df.groupby(['Genus'])[comp_class].transform('mean')
    genera_df['expected_total_mean'] = expected_mean

    genera_df = genera_df[['Genus', 'tested_compounds_count', 'hit_compounds_count', 'mean_compound_value', 'expected_total_mean']]
    genera_df = genera_df.reset_index(drop=True)

    genera_df = genera_df.drop_duplicates(subset=['Genus'])
    genera_df.reset_index(drop=True, inplace=True)
    assert len(genera_df) == num_genera_tested
    genera_df.to_csv(os.path.join(processed_compounds_output_path, f'genus_level_{comp_class}_information.csv'))
    return genera_df


def get_chemical_diversity_measure_for_genera():
    pass


def main():
    # get_processed_metabolite_data()
    add_class_information()
    for comp_class in CLASSES_OF_INTEREST:
        get_genus_level_version(comp_class)


if __name__ == '__main__':
    main()
