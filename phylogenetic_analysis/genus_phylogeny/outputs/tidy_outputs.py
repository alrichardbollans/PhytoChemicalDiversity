import os

import pandas as pd
from phytochempy.compound_properties import NP_PATHWAYS


def holm_correction(df: pd.DataFrame, p_value_col: str):
    # holm_simple_1979
    # Sture Holm, ‘A Simple Sequentially Rejective Multiple Test Procedure’, Scandinavian Journal of Statistics 6, no. 2 (1979): 65–70.

    # The adjustment is the same as hochberg, but the usage is slightly different.
    new_df = df.sort_values(by=p_value_col)
    n = len(new_df.index)
    new_df.reset_index(inplace=True, drop=True)

    new_df['holm_adjusted_p_value'] = new_df.apply(lambda x: x[p_value_col] * (n - x.name), axis=1)

    return new_df


def diversity():
    indices = ['FAD', 'MFAD', 'APWD', 'Hbc', 'H', 'G']
    output = []
    for ind in indices:
        result = pd.read_csv(os.path.join(ind, 'Genus_phylogenetic_signal_results.csv'), index_col=0)
        result['index'] = ind
        result = result[['index'] + [c for c in result.columns if c != 'index']]
        output.append(result)
    final_output = pd.concat(output)
    holm_correction(final_output, 'pvalue').to_csv('diversity_results.csv')

# def exploration():
#     indices = ['exploration_index', 'species_richness', 'N']
#     output = []
#     for ind in indices:
#         result = pd.read_csv(os.path.join(ind, 'Genus_phylogenetic_signal_results.csv'), index_col=0)
#         result['index'] = ind
#         result = result[['index'] + [c for c in result.columns if c != 'index']]
#         output.append(result)
#     final_output = pd.concat(output)
#     holm_correction(final_output, 'pvalue').to_csv('exploration_results.csv')

def pathways():
    output1 = []
    output2 = []
    for p in NP_PATHWAYS:
        result1 = pd.read_csv(os.path.join(f'mean_identified_as_{p.replace(" ","_")}', 'Genus_phylogenetic_signal_results.csv'), index_col=0)
        result2 = pd.read_csv(os.path.join(f'norm_mean_identified_as_{p.replace(" ","_")}', 'Genus_phylogenetic_signal_results.csv'), index_col=0)
        result1['pathway'] = p.replace(' ','_')
        result2['pathway'] = p.replace(' ','_')
        result1 = result1[['pathway'] + [c for c in result1.columns if c != 'index']]
        result2 = result2[['pathway'] + [c for c in result2.columns if c != 'index']]
        output1.append(result1)
        output2.append(result2)

    holm_correction(pd.concat(output1), 'pvalue').to_csv('pathway_results.csv')
    holm_correction(pd.concat(output2), 'pvalue').to_csv('norm_pathway_results.csv')


def main():
    diversity()
    pathways()


if __name__ == '__main__':
    main()
