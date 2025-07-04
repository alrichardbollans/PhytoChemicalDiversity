import os

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from phytochempy.compound_properties import get_npclassifier_pathway_columns_in_df, NP_PATHWAYS

from collect_and_compile_data.collect_compound_data import all_species_compound_csv, COMPOUND_ID_COL

out_dir = os.path.join('outputs', 'plots')


def histogram_of_count_per_species():
    data_df = compound_data[['accepted_species', COMPOUND_ID_COL]].drop_duplicates()
    group_size = data_df.groupby(['accepted_species']).size()
    pd.DataFrame(group_size).to_csv(os.path.join(out_dir, 'counts_per_species.csv'))
    group_size = group_size[group_size <= 100]
    sns.histplot(x=group_size, kde=True)
    plt.xlabel('Number of Compounds')
    plt.savefig(os.path.join(out_dir, 'histogram_of_count_per_species_lt_100.png'))
    plt.close()

    group_size = data_df.groupby(['accepted_species']).size()
    group_size = group_size[group_size > 100]
    sns.histplot(x=group_size, kde=True)
    plt.xlabel('Number of Compounds')
    plt.savefig(os.path.join(out_dir, 'histogram_of_count_per_species_gt_100.png'))
    plt.close()


def histogram_of_count_per_pathway():
    sns.set_theme(style='white')

    data_df = compound_data
    pathway_cols = get_npclassifier_pathway_columns_in_df(data_df)
    data_df = compound_data[[COMPOUND_ID_COL] + pathway_cols].drop_duplicates(subset=COMPOUND_ID_COL)
    pathway_counts = {}
    for pathway in NP_PATHWAYS:
        pathway_counts[pathway] = 0
        for p in pathway_cols:
            pathway_data_count = len(data_df[data_df[p] == pathway])
            pathway_counts[pathway] += pathway_data_count

    counts = pd.DataFrame(pathway_counts, index=['Number of Compounds']).T
    counts.columns = ['Number of Compounds']
    counts['Pathway'] = counts.index
    counts = counts.sort_values(by='Number of Compounds', ascending=False)
    counts.to_csv(os.path.join(out_dir, 'counts_per_pathway.csv'))
    counts['Pathway'] = ['a.', 'b.', 'c.', 'd.', 'e.', 'f.', 'g.']
    # plotting data on chart
    print(counts)
    sns.barplot(counts.reset_index(), x="Pathway", y="Number of Compounds", hue='Pathway')
    plt.legend()
    plt.savefig(os.path.join(out_dir, 'counts_for_pathways.png'))
    plt.close()


def main():
    histogram_of_count_per_species()
    histogram_of_count_per_pathway()


if __name__ == '__main__':
    compound_data = pd.read_csv(all_species_compound_csv)
    main()
