import os

import pandas as pd
from pkg_resources import resource_filename

from collect_and_compile_data.get_diversity_metrics.gather_diversity_measures import genus_diversity_data_csv

_output_path = resource_filename(__name__, 'outputs')


def main():

    diversity_df = pd.read_csv(genus_diversity_data_csv, index_col=0)
    indices = ['FAD','MFAD', 'APWD', 'H', 'Hbc', 'G', 'J']
    # Top 5 for each index
    for i in indices:
        top_5 = diversity_df.sort_values(by=i, ascending=False).head(5)
        top_5.to_csv(os.path.join(_output_path, 'top_5_' + i + '.csv'))

    # Scale to make means comparable
    from sklearn.preprocessing import StandardScaler
    for index in indices:
        scaler = StandardScaler()
        diversity_df[index + '_scaled'] = scaler.fit_transform(diversity_df[[index]])

    diversity_df['pathway_mean'] = diversity_df[['H_scaled', 'Hbc_scaled', 'G_scaled', 'J_scaled']].mean(axis=1)
    top_5_mean = diversity_df.sort_values(by='pathway_mean', ascending=False).head(5)
    top_5_mean.to_csv(os.path.join(_output_path, 'top_5_pathway_mean.csv'))

    diversity_df['distance_mean'] = diversity_df[['FAD_scaled','MFAD_scaled', 'APWD_scaled']].mean(axis=1)
    top_5_mean = diversity_df.sort_values(by='distance_mean', ascending=False).head(5)
    top_5_mean.to_csv(os.path.join(_output_path, 'top_5_distance_mean.csv'))

    # Top 5 mean of indices
    diversity_df['mean'] = diversity_df[
        ['H_scaled', 'Hbc_scaled', 'G_scaled', 'J_scaled', 'FAD_scaled','MFAD_scaled', 'APWD_scaled']].mean(axis=1)
    top_5_mean = diversity_df.sort_values(by='mean', ascending=False).head(5)
    top_5_mean.to_csv(os.path.join(_output_path, 'top_5_mean.csv'))


if __name__ == '__main__':
    main()
