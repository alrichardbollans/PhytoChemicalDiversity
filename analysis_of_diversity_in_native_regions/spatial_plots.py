import os
from pathlib import Path

from pkg_resources import resource_filename

from analysis_of_diversity_in_native_regions.analyse_relation_to_pd import get_working_data
from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import FAD_INDICES, PATHWAY_INDICES, METRICS

_inputs_path = resource_filename(__name__, 'inputs')

_temp_outputs_path = resource_filename(__name__, 'temp_outputs')

_output_path = resource_filename(__name__, 'outputs')


def plot_dist_of_metric(df_with_region_data, metric, colormap: str = 'viridis', out_path: str = None):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import cartopy.io.shapereader as shpreader

    tdwg3_shp = shpreader.Reader(
        os.path.join(_inputs_path, 'wgsrpd-master', 'level3', 'level3.shp'))
    tdwg3_region_codes = df_with_region_data['Group'].values

    ## Colour maps range is 0 - 1, so the values are standardised for this
    max_val = df_with_region_data[metric].max()
    min_val = df_with_region_data[metric].min()
    norm = plt.Normalize(min_val, max_val)
    print('plotting countries')

    plt.figure(figsize=(15, 9.375))
    ax = plt.axes(projection=ccrs.Mollweide())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, linewidth=2)

    cmap = mpl.colormaps[colormap]
    for country in tdwg3_shp.records():

        tdwg_code = country.attributes['LEVEL3_COD']
        if tdwg_code in tdwg3_region_codes:
            ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                              facecolor=cmap(
                                  norm(df_with_region_data.loc[df_with_region_data['Group'] == tdwg_code, metric].iloc[
                                           0])),
                              label=tdwg_code)

        else:
            ax.add_geometries([country.geometry], ccrs.PlateCarree(),
                              facecolor='white',
                              label=tdwg_code)

    all_map_isos = [country.attributes['LEVEL3_COD'] for country in tdwg3_shp.records()]
    missed_names = [x for x in tdwg3_region_codes if x not in all_map_isos]
    print(f'iso codes not plotted on map: {missed_names}')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    plt.tight_layout()
    fig = plt.gcf()
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.175, 0.02, 0.65])
    cbar1 = fig.colorbar(sm, cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=30)

    if out_path is None:
        out_dir = os.path.join('outputs', 'spatial_plots')
        if metric in FAD_INDICES:
            out_dir = os.path.join(out_dir, 'functional_indices')
        if metric in PATHWAY_INDICES:
            out_dir = os.path.join(out_dir, 'pathway_indices')

        out_path = os.path.join(out_dir, f'{metric}.jpg')
    Path(os.path.dirname(out_path)).mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=400, bbox_inches='tight')
    plt.close()
    plt.cla()
    plt.clf()


def main():
    working_data = get_working_data()
    plot_dist_of_metric(working_data, 'number_of_species_in_group')
    plot_dist_of_metric(working_data, 'phylogenetic_diversity')

    for metric in METRICS:
        plot_dist_of_metric(working_data, metric)


if __name__ == '__main__':
    main()
