import os

import pandas as pd
from pkg_resources import resource_filename
from rdkit.Chem import PandasTools
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.DataManip.Metric import GetTanimotoDistMat

from library_info_and_data_import import get_processed_metabolite_data

_inputs_path = resource_filename(__name__, 'inputs')

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')

_output_path = resource_filename(__name__, 'outputs')
if not os.path.isdir(_inputs_path):
    os.mkdir(_inputs_path)
if not os.path.isdir(_temp_outputs_path):
    os.mkdir(_temp_outputs_path)
if not os.path.isdir(_output_path):
    os.mkdir(_output_path)


def read_data(taxon_grouping: str):
    all_taxa_metabolite_data = get_processed_metabolite_data('SMILES', taxon_grouping=taxon_grouping)[['Genus', 'SMILES']]

    return all_taxa_metabolite_data


def get_pairwise_distances_from_data(df: pd.DataFrame):
    PandasTools.AddMoleculeColumnToFrame(df, 'SMILES', 'Molecule', includeFingerprints=True)
    df = df.dropna(subset=['Molecule'])[['Molecule', 'SMILES']]

    # Produce a hashed Morgan fingerprint for each molecule
    df['morgan_fingerprint'] = df['Molecule'].apply(lambda x: GetMorganFingerprintAsBitVect(x, 2))
    distmat = GetTanimotoDistMat(df['morgan_fingerprint'].values)

    return distmat


def calculate_FAD_measures(df: pd.DataFrame, taxon_grouping: str = 'Genus'):
    FAD_outputs = {}
    MFAD_outputs = {}
    APWD_outputs = {}
    N_outputs = {}
    for taxon in df[taxon_grouping].unique():
        taxon_data = df[df[taxon_grouping] == taxon]
        if len(taxon_data) > 1:
            distances = get_pairwise_distances_from_data(taxon_data)
            FAD_outputs[taxon] = distances.sum()
            MFAD_outputs[taxon] = distances.sum() / len(taxon_data)
            APWD_outputs[taxon] = distances.sum() / len(distances)
            N_outputs[taxon] = len(taxon_data)

    out_df = pd.DataFrame.from_dict(FAD_outputs, orient='index', columns=['FAD'])
    out_df['MFAD'] = MFAD_outputs.values()
    out_df['APWD'] = APWD_outputs.values()
    out_df['N'] = N_outputs.values()
    out_df = out_df.reset_index(names=[taxon_grouping])
    return out_df


def main():
    metabolite_data = read_data('Genus')
    FAD_measures = calculate_FAD_measures(metabolite_data, taxon_grouping='Genus')
    FAD_measures.to_csv(os.path.join('outputs','genus_level_distance_diversity_information.csv'))

if __name__ == '__main__':
    main()
