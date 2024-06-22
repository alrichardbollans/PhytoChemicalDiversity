import os

import pandas as pd
from pkg_resources import resource_filename

from collect_compound_data import FAMILIES_OF_INTEREST

_output_path = resource_filename(__name__, 'outputs')
species_richness_csv = os.path.join(_output_path, 'species_richness.csv')
if not os.path.isdir(_output_path):
    os.mkdir(_output_path)

if __name__ == '__main__':
    from wcvpy.wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns

    all_taxa = get_all_taxa(families_of_interest=FAMILIES_OF_INTEREST, accepted=True)

    species_df = all_taxa[all_taxa[wcvp_columns['rank']] == 'Species']
    counted = species_df.groupby([wcvp_accepted_columns['parent_name']]).size()
    counted_df = pd.DataFrame({'Genus': counted.index, 'species_richness': counted.values})
    counted_df.to_csv(species_richness_csv)
