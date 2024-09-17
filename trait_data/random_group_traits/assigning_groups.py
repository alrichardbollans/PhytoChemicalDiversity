import os
from random import choice

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import hist


def main():
    species_data = pd.read_csv(os.path.join('..', 'compile_trait_data', 'outputs', 'species_trait_data.csv'))[
        ['accepted_species', 'Genus','Animal Richness', 'PC0', 'PC1']]

    # For the moment, just mimic number of genera.
    number_of_groups = len(species_data['Genus'].unique().tolist())

    groups = range(0, number_of_groups)
    assigned_groups = [str(choice(groups)) for i in species_data.index]
    # hist(assigned_groups)
    # plt.show()
    species_data['Assigned_group'] = assigned_groups

    # After groups have been assigned, add compound data
    compound_data = pd.read_csv(os.path.join('..', 'collect_compound_data', 'outputs', 'all_species_compound_data.csv'))
    working_data = pd.merge(compound_data, species_data, how='left', on='accepted_species')
    # then need to calculate metrics.



    pass


if __name__ == '__main__':
    main()
