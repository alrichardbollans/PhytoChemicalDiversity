import os

import pandas as pd


from library_info_and_data_import import processed_compounds_output_path




def your_function():
    diversity_df = pd.read_csv(os.path.join(processed_compounds_output_path, f'genus_level_pathway_diversity_information.csv'))
    indices = ['shannon_index','gini_index','normalised_gini','normalised_shannon']
    for i in indices:
        for j in indices:
            if i != j:
                correlation = diversity_df[[i, j]].corr().iloc[0, 1]
                print(f"Correlation between {i} and {j}: {correlation}")



def main():
    your_function()


if __name__ == '__main__':
    main()
