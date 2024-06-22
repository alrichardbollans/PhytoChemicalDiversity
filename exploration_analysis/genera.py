import io
import os

import folium
import pandas as pd
import seaborn
from PIL import Image
from matplotlib import pyplot as plt
from pkg_resources import resource_filename
from wcvpy.wcvp_download import plot_native_number_accepted_taxa_in_regions, get_distributions_for_accepted_taxa
from wcvpy.wcvp_download.plot_distributions import _OHE_native_dists

from collect_compound_data import genus_pathway_data_csv
from diversity_metrics import species_richness_csv

_map_params = {'location': [40, 20], 'zoom_start': 2, 'font_size': '1.5rem', 'zoom_control': False, 'tiles':None}

_output_path = resource_filename(__name__, 'outputs')


def setup():
    genus_df = pd.read_csv(genus_pathway_data_csv, index_col=0)
    genus_df = genus_df.rename(columns={'identified_compounds_count': 'N'})
    richness_df = pd.read_csv(species_richness_csv, index_col=0)
    sampling_df = pd.merge(richness_df, genus_df, on='Genus', how='left')[['Genus', 'species_richness', 'N']]
    sampling_df['N'] = sampling_df['N'].fillna(0)
    sampling_df['exploration_index'] = sampling_df['N'] / sampling_df['species_richness']
    sampling_df['norm_exploration_index'] = sampling_df['exploration_index'] / sampling_df['exploration_index'].max()
    sampling_df['norm_N'] = sampling_df['N'] / sampling_df['N'].max()
    sampling_df['norm_species_richness'] = sampling_df['species_richness'] / sampling_df['species_richness'].max()
    sampling_df.to_csv(os.path.join(_output_path, 'sampling_effort.csv'))
    return sampling_df


def sampling_effort_plot():
    sampling_df = setup()
    seaborn.regplot(data=sampling_df, y='N', x='species_richness')

    plt.xlabel('Species Richness')
    plt.ylabel('Number of Identified Compounds')
    plt.savefig(os.path.join('outputs', 'sampling_effort.jpg'), dpi=300)
    plt.close()


def exploration_index_comparison():
    sampling_df = setup()

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


def get_global_distribution_of_exploration_index_and_N():
    sampling_df = setup()
    df_with_dists = get_distributions_for_accepted_taxa(sampling_df.drop_duplicates(subset=['Genus'], keep='first'), 'Genus')
    ohe = _OHE_native_dists(df_with_dists)
    # For each region, look at the genera in that region and work out the mean exploration index
    # Then plot on map
    regions = [c for c in ohe.columns if c not in df_with_dists.columns]
    out_dict = {}
    for region in regions:
        region_df = ohe[ohe[region] == 1][['Genus', 'exploration_index', region]]
        out_dict[region] = [region_df['exploration_index'].mean()]

    out_df = pd.DataFrame.from_dict(out_dict, columns=['mean_exploration_index'], orient='index').reset_index(names=['Region'])
    out_df.to_csv(os.path.join('outputs', 'global_distribution_of_exploration_index.csv'))

    # For each region, look at the genera in that region and work out the mean
    # Then plot on map

    out_dict = {}
    for region in regions:
        region_df = ohe[ohe[region] == 1][['Genus', 'N', region]]
        out_dict[region] = [region_df['N'].mean()]

    out_df = pd.DataFrame.from_dict(out_dict, columns=['mean_N'], orient='index').reset_index(names=['Region'])
    out_df.to_csv(os.path.join('outputs', 'global_distribution_of_N.csv'))

    # For each region, look at the genera in that region and work out the mean
    # Then plot on map

    out_dict = {}
    for region in regions:
        region_df = ohe[ohe[region] == 1][['Genus', 'species_richness', region]]
        out_dict[region] = [region_df['species_richness'].mean()]

    out_df = pd.DataFrame.from_dict(out_dict, columns=['mean_species_richness'], orient='index').reset_index(names=['Region'])
    out_df.to_csv(os.path.join('outputs', 'global_distribution_of_species_richness_for_genera.csv'))


def plot_global_distribution_of_exploration_index():
    import geopandas as gpd
    index_df = pd.read_csv(os.path.join('outputs', 'global_distribution_of_exploration_index.csv'), index_col=0)
    # Now plot
    # Create the choropleth map
    m = folium.Map(**_map_params)
    world = gpd.read_file('inputs/wgsrpd-master/level3/level3.shp')
    folium.Choropleth(
        geo_data=world,
        use_jenks=True,
        name='Mean Genus Exploration Index',
        data=index_df,
        columns=['Region', 'mean_exploration_index'],
        key_on='feature.properties.LEVEL3_COD',
        fill_color='BuPu', highlight=True,
        fill_opacity=0.8,
        line_opacity=0.2,
        overlay=False,
        legend_name='Mean Exploration Index'
    ).add_to(m)

    img_data = m._to_png(2)
    img = Image.open(io.BytesIO(img_data))
    img.save('outputs/genus_exploration_index_map.png')

    # Display the map
    m.save('outputs/genus_exploration_index_map.html')


def plot_global_distribution_of_N():
    import geopandas as gpd
    index_df = pd.read_csv(os.path.join('outputs', 'global_distribution_of_N.csv'), index_col=0)
    # Now plot
    # Create the choropleth map
    m = folium.Map(**_map_params)
    world = gpd.read_file('inputs/wgsrpd-master/level3/level3.shp')
    folium.Choropleth(
        geo_data=world,
        use_jenks=True,
        name='Mean Number of Identified Compounds',
        data=index_df,
        columns=['Region', 'mean_N'],
        key_on='feature.properties.LEVEL3_COD',
        fill_color='BuPu', highlight=True,
        fill_opacity=0.8,
        line_opacity=0.2,
        overlay=False,
        legend_name='Mean Number of Identified Compounds'
    ).add_to(m)

    img_data = m._to_png(2)
    img = Image.open(io.BytesIO(img_data))
    img.save('outputs/genus_N_map.png')

    # Display the map
    m.save('outputs/genus_N_map.html')


def plot_global_distribution_of_sp_richness():
    import geopandas as gpd
    index_df = pd.read_csv(os.path.join('outputs', 'global_distribution_of_species_richness_for_genera.csv'), index_col=0)
    # Now plot
    # Create the choropleth map
    m = folium.Map(**_map_params)
    world = gpd.read_file('inputs/wgsrpd-master/level3/level3.shp')
    folium.Choropleth(
        geo_data=world,
        use_jenks=True,
        name='Mean Number of Species',
        data=index_df,
        columns=['Region', 'mean_species_richness'],
        key_on='feature.properties.LEVEL3_COD',
        fill_color='BuPu', highlight=True,
        fill_opacity=0.8,
        line_opacity=0.2,
        overlay=False,
        legend_name='Mean Number of Species'
    ).add_to(m)

    img_data = m._to_png(2)
    img = Image.open(io.BytesIO(img_data))
    img.save('outputs/genus_species_richness_map.png')

    # Display the map
    m.save('outputs/genus_species_richness_map.html')


def main():
    sampling_effort_plot()
    exploration_index_comparison()
    get_global_distribution_of_exploration_index_and_N()
    plot_global_distribution_of_exploration_index()
    plot_global_distribution_of_N()
    plot_global_distribution_of_sp_richness()


if __name__ == '__main__':
    main()
