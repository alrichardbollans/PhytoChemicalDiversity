import math
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
# Uniqueness is up to 'InChIKey'
# TODO: Check usage of full key when collecting data
COMPOUND_ID_COL = 'InChIKey'
FAMILIES_OF_INTEREST = ['Apocynaceae', 'Loganiaceae', 'Rubiaceae']
CLASSES_OF_INTEREST = ['alkaloid', 'steroid', 'terpenoid', 'peptide', 'lipid', 'polyketide']
MINOR_CLASSES = ['pyrrolizidine', 'quinoline']
output_compound_info = [wcvp_accepted_columns['species_w_author'], 'Genus',
                        'example_metabolite_name', 'InChIKey', 'SMILES', 'InChIKey_simp', 'active_chembl_compound',
                        'mean_ic50_μM', 'lipinski_pass', 'veber_pass',
                        'MAIP_model_score',
                        ] + compound_class_columns + [wcvp_accepted_columns['species'], wcvp_accepted_columns['family']]


def get_processed_metabolite_data():
    ## First get assay data
    chembl_compound_assays = pd.read_csv(os.path.join(chembl_apm_assay_info_csv),
                                         index_col=0).dropna(subset=[COMPOUND_ID_COL, 'assay_pchembl_value'], how='any')

    chembl_compound_assays_summary = chembl_compound_assays[[COMPOUND_ID_COL, 'mean_ic50_μM']].drop_duplicates(keep='first')
    assert len(chembl_compound_assays_summary[COMPOUND_ID_COL].unique().tolist()) == len(chembl_compound_assays[COMPOUND_ID_COL].unique().tolist())

    chembl_compound_assays_summary = chembl_compound_assays_summary.sort_values(by=COMPOUND_ID_COL).reset_index(drop=True)
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

    class_suggestions = processed_metabolite_data['NPclassif_pathway_results'].unique()
    print('NPclassif_pathway_results class suggestions')
    print(class_suggestions)

    class_suggestions = processed_metabolite_data['classyfire_superclass'].unique()
    print('classyfire_superclass class suggestions')
    print(class_suggestions)


def add_class_information():
    metabolite_data = pd.read_csv(_processed_metabolite_csv, index_col=0)
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
    metabolite_data.describe(include='all').to_csv(os.path.join(processed_compounds_output_path, 'processed_metabolite_summary.csv'))


def get_genus_level_version_for_compound(comp_class: str):
    metabolite_data = pd.read_csv(processed_metabolite_with_classes_csv, index_col=0)
    expected_mean = metabolite_data[comp_class].mean()

    genera_df = metabolite_data.copy()

    num_genera_tested = len(genera_df['Genus'].unique().tolist())

    # count labelled species
    counts = genera_df.value_counts('Genus')
    counts.name = 'identified_compounds_count'
    genera_df = pd.merge(genera_df, counts, how='left', left_on='Genus', right_index=True)
    genera_df['identified_comp_class_count'] = genera_df.groupby(['Genus'])[comp_class].transform('sum')
    genera_df['mean_identified_as_class'] = genera_df.groupby(['Genus'])[comp_class].transform('mean')
    genera_df['expected_total_mean'] = expected_mean

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

    genera_df['weigthing_factor'] = genera_df['identified_compounds_count'].apply(weighting_factor)

    genera_df['norm_mean_identified_as_class'] = (genera_df['weigthing_factor'] * genera_df['mean_identified_as_class']) + (
            1 - genera_df['weigthing_factor']) * expected_mean

    genera_df = genera_df[
        ['Genus', 'identified_compounds_count', 'identified_comp_class_count', 'mean_identified_as_class', 'expected_total_mean', 'weigthing_factor',
         'norm_mean_identified_as_class']]
    genera_df = genera_df.reset_index(drop=True)

    genera_df = genera_df.drop_duplicates(subset=['Genus'])
    genera_df.reset_index(drop=True, inplace=True)
    assert len(genera_df) == num_genera_tested
    genera_df.to_csv(os.path.join(processed_compounds_output_path, f'genus_level_{comp_class}_information.csv'))
    return genera_df


def get_chemical_diversity_measure_for_genera():
    metabolite_data = pd.read_csv(processed_metabolite_with_classes_csv, index_col=0)
    genera_df = metabolite_data[['Genus']]
    # count labelled species
    counts = genera_df.value_counts('Genus')
    counts.name = 'identified_compounds_count'
    genera_df = pd.merge(genera_df, counts, how='left', left_on='Genus', right_index=True)
    genera_df = genera_df.drop_duplicates(subset=['Genus'])
    num_genera_tested = len(genera_df['Genus'].unique().tolist())

    for comp_class in CLASSES_OF_INTEREST:
        genera_comp_df = metabolite_data.copy()

        genera_comp_df[comp_class] = genera_comp_df.groupby(['Genus'])[comp_class].transform('max')

        genera_comp_df = genera_comp_df[['Genus', comp_class]]
        genera_comp_df = genera_comp_df.reset_index(drop=True)

        genera_comp_df = genera_comp_df.drop_duplicates(subset=['Genus'])
        genera_comp_df.reset_index(drop=True, inplace=True)
        assert len(genera_comp_df) == num_genera_tested
        genera_df = pd.merge(genera_df, genera_comp_df, how='left', on='Genus')
        assert len(genera_df) == num_genera_tested
    ### This is just a crude measure for thhe moment
    genera_df['norm_factor'] = genera_df['identified_compounds_count'].apply(lambda x: min(x, len(CLASSES_OF_INTEREST)))
    genera_df['mean_identified_as_class'] = genera_df[CLASSES_OF_INTEREST].sum(axis=1)
    genera_df['mean_identified_as_class'] = genera_df['mean_identified_as_class'] / genera_df['norm_factor']
    genera_df['mean_identified_as_class'] = genera_df['mean_identified_as_class'].apply(lambda x: min(x, 1))

    genera_df.to_csv(os.path.join(processed_compounds_output_path, f'genus_level_chemical_diversity_information.csv'))


def main():
    # get_processed_metabolite_data()
    # add_class_information()
    for comp_class in CLASSES_OF_INTEREST + MINOR_CLASSES:
        get_genus_level_version_for_compound(comp_class)

    get_chemical_diversity_measure_for_genera()


if __name__ == '__main__':
    main()
