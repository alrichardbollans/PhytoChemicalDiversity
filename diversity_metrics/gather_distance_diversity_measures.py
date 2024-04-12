import os

import pandas as pd
from pkg_resources import resource_filename
from rdkit.Chem import PandasTools
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.DataManip.Metric import GetTanimotoDistMat

from collect_compound_data import get_relevant_deduplicated_data, COMPOUND_ID_COL, all_taxa_compound_csv, FAMILIES_OF_INTEREST

_output_path = resource_filename(__name__, 'outputs')
genus_distance_diversity_data_csv = os.path.join(_output_path, 'genus_level_distance_diversity_information.csv')
if not os.path.isdir(_output_path):
    os.mkdir(_output_path)


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
        else:
            FAD_outputs[taxon] = 0
            MFAD_outputs[taxon] = 0
            APWD_outputs[taxon] = 0
            N_outputs[taxon] = len(taxon_data)

    out_df = pd.DataFrame.from_dict(FAD_outputs, orient='index', columns=['FAD'])
    out_df['norm_FAD'] = out_df['FAD'] / out_df['FAD'].max()
    out_df['MFAD'] = MFAD_outputs.values()
    out_df['norm_MFAD'] = out_df['MFAD']/out_df['MFAD'].max()
    out_df['APWD'] = APWD_outputs.values()
    out_df['N'] = N_outputs.values()
    out_df = out_df.reset_index(names=[taxon_grouping])
    return out_df


def main():
    my_df = pd.read_csv(all_taxa_compound_csv, index_col=0)
    processed = get_relevant_deduplicated_data(my_df, COMPOUND_ID_COL, 'Genus', FAMILIES_OF_INTEREST)
    FAD_measures = calculate_FAD_measures(processed, taxon_grouping='Genus')
    FAD_measures.to_csv(genus_distance_diversity_data_csv)


if __name__ == '__main__':
    main()
