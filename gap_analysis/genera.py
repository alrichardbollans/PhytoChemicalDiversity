import os

import pandas as pd
import seaborn
from matplotlib import pyplot as plt
from pkg_resources import resource_filename
from wcvp_download import plot_native_number_accepted_taxa_in_regions
from wcvp_name_matching import get_accepted_info_from_names_in_column

from collect_compound_data import genus_pathway_data_csv
from diversity_metrics import genus_distance_diversity_data_csv, species_richness_csv

_output_path = resource_filename(__name__, 'outputs')


def setup():
    genus_df = pd.read_csv(genus_pathway_data_csv, index_col=0)
    genus_df = genus_df.rename(columns={'identified_compounds_count': 'N'})
    richness_df = pd.read_csv(species_richness_csv, index_col=0)
    sampling_df = pd.merge(richness_df, genus_df, on='Genus', how='left')[['Genus', 'species_richness', 'N']]
    sampling_df['N'] = sampling_df['N'].fillna(0)
    return sampling_df


def sampling_effort_plot():
    sampling_df = setup()
    seaborn.regplot(data=sampling_df, y='N', x='species_richness')

    plt.xlabel('Species Richness')
    plt.ylabel('Number of Identified Compounds')
    plt.savefig(os.path.join('outputs', 'sampling_effort.jpg'), dpi=300)
    plt.close()


def exploration_index():
    sampling_df = setup()
    sampling_df['exploration_index'] = sampling_df['N'] / sampling_df['species_richness']
    # most_sampled = sampling_df.sort_values(by='species_richness', ascending=False)
    most_sampled = sampling_df.sort_values(by=['exploration_index', 'species_richness'], ascending=False)
    most_sampled.to_csv(os.path.join('outputs', 'most_sampled.csv'))

    plot_native_number_accepted_taxa_in_regions(most_sampled.head(20), 'Genus', _output_path,
                                                'most_sampled.jpg', include_extinct=False)

    # Plot the top 50 for most and least sampled.
    least_sampled = sampling_df.sort_values(by=['exploration_index', 'species_richness'], ascending=[True, False])
    least_sampled.to_csv(os.path.join('outputs', 'least_sampled.csv'))
    plot_native_number_accepted_taxa_in_regions(least_sampled.head(20), 'Genus', _output_path,
                                                'least_sampled.jpg', include_extinct=False)

    pass


def main():
    sampling_effort_plot()
    exploration_index()


if __name__ == '__main__':
    main()
