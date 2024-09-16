import os

import pandas as pd

from get_phylogeny.genus_phylogeny.outputs.tidy_outputs import holm_correction


def correlations():
    tax_div = pd.read_csv('correlations_with_num_species.csv', index_col=0)
    tax_div['var'] = 'Taxonomic Diversity'
    phy_div = pd.read_csv('correlations_with_phyl_diversity.csv', index_col=0)
    phy_div['var'] = 'Phylogenetic Diversity'

    final_output = pd.concat([tax_div, phy_div])
    holm_correction(final_output, 'pvalue').to_csv('correlation_results.csv')


def pgls():
    tax_div = pd.read_csv(os.path.join('pgls', 'p_values.csv'), index_col=0)
    transposed = tax_div.T
    transposed['metric'] = transposed.index
    tax_div = transposed[['metric', 'Num. Species']]
    phy_div = transposed[['metric', 'Phylo diversity']]

    tax_div = tax_div.rename(columns={'Num. Species': 'pvalue'})
    phy_div = phy_div.rename(columns={'Phylo diversity': 'pvalue'})

    tax_div['var'] = 'Taxonomic Diversity'
    phy_div['var'] = 'Phylogenetic Diversity'

    final_output = pd.concat([tax_div, phy_div])
    holm_correction(final_output, 'pvalue').to_csv(os.path.join('pgls', 'corrected_p_values.csv'))


def main():
    correlations()
    pgls()


if __name__ == '__main__':
    main()
